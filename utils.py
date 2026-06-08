import asyncio
import itertools
import os
import re
import shutil
from asyncio.subprocess import Process
from dataclasses import dataclass
from enum import Enum
from pathlib import Path
from tempfile import mkdtemp
from typing import Protocol

from Bio import SeqIO
from pyfixtures import fixture
from structlog import get_logger
from virtool.caches.utils import derive_key
from virtool.models.enums import LibraryType
from virtool.workflow import RunSubprocess
from virtool.workflow.analysis import ReadPaths
from virtool.workflow.data.cache import CacheHit, WorkflowCache
from virtool.workflow.data.samples import WFSample
from virtool.workflow.utils import get_workflow_version

BOWTIE2_BUILD_TOOL = "bowtie2-build"
SKEWER_TOOL = "skewer"
WORKFLOW_NAME = "nuvs"
WORKFLOW_VERSION = get_workflow_version()

logger = get_logger("workflow")


class SkewerMode(str, Enum):
    """The mode to run Skewer in."""

    PAIRED_END = "pe"
    """Run Skewer in paired-end mode."""

    SINGLE_END = "any"
    """Run Skewer in single-end mode."""


@dataclass
class SkewerConfiguration:
    """A configuration for running Skewer."""

    min_length: int
    """The minimum length of a trimmed read."""

    mode: SkewerMode
    """The mode to run Skewer in."""

    end_quality: int = 20
    """The minimum quality score for the end of a trimmed read."""

    max_error_rate: float = 0.1
    """
    The maximum error rate for a trimmed read. Reads that exceed the rate will be
    discarded.
    """

    max_indel_rate: float = 0.03
    """
    The maximum indel rate for a trimmed read. Reads that exceed the rate will be
    discarded.
    """

    mean_quality: int = 25
    """The minimum mean quality score for a trimmed read. Reads  """

    number_of_processes: int = 1
    """The number of processes to use when running Skewer."""

    quiet: bool = True
    """Whether to run Skewer in quiet mode."""

    other_options: tuple[str] = ("-n", "-z")
    """Other options to pass to Skewer."""


@dataclass
class SkewerResult:
    """Represents the result of running Skewer to trim paired or unpaired FASTQ data."""

    command: list[str]
    """The command used to run Skewer."""

    output_path: Path
    """The path to the directory containing the trimmed reads."""

    process: Process
    """The process that ran Skewer."""

    read_paths: ReadPaths
    """The paths to the trimmed reads."""

    @property
    def left(self) -> Path:
        """The path to the left reads of paired or unpaired FASTQ data."""
        return self.read_paths[0]

    @property
    def right(self) -> Path | None:
        """The path to the right reads of a paired Illumina dataset.

        Set to ``None`` if the dataset in unpaired.
        """
        try:
            return self.read_paths[1]
        except IndexError:
            return None


class SkewerRunner(Protocol):
    """A protocol describing callables that can be used to run Skewer."""

    async def __call__(
        self,
        config: SkewerConfiguration,
        paths: ReadPaths,
        output_path: Path,
    ) -> SkewerResult: ...


def calculate_skewer_trimming_parameters(
    sample: WFSample,
    min_read_length: int,
) -> SkewerConfiguration:
    """Calculates trimming parameters based on the library type, and minimum allowed trim length.

    :param sample: The sample to calculate trimming parameters for.
    :param min_read_length: The minimum length of a read before it is discarded.
    :return: the trimming parameters
    """
    config = SkewerConfiguration(
        min_length=min_read_length,
        mode=SkewerMode.PAIRED_END if sample.paired else SkewerMode.SINGLE_END,
    )

    if sample.library_type == LibraryType.amplicon:
        config.end_quality = 0
        config.mean_quality = 0
        config.min_length = min_read_length

        return config

    if sample.library_type == LibraryType.srna:
        config.max_length = 22
        config.min_length = 20

        return config

    raise ValueError(f"Unknown library type: {sample.library_type}")


@fixture
def skewer(proc: int, run_subprocess: RunSubprocess) -> SkewerRunner:
    """Provides an asynchronous function that can run skewer.

    The provided function takes a :class:`.SkewerConfiguration` and a tuple of paths to
    the left and right reads to trim. If a single member tuple is provided, the dataset
    is assumed to be unpaired.

    The Skewer process will automatically be assigned the number of processes configured
    for the workflow run.

    Example:
    -------
    .. code-block:: python

        @step
        async def step_one(skewer: SkewerRunner, work_path: Path):
           config = SkewerConfiguration(
               mean_quality=30
           )

           skewer_result = await skewer(config, (
               work_path / "test_1.fq.gz",
               work_path / "test_2.fq.gz",
           ))


    """
    if shutil.which("skewer") is None:
        raise RuntimeError("skewer is not installed.")

    async def func(
        config: SkewerConfiguration,
        read_paths: ReadPaths,
        output_path: Path,
    ):
        temp_path = Path(await asyncio.to_thread(mkdtemp, suffix="_virtool_skewer"))

        await asyncio.to_thread(output_path.mkdir, exist_ok=True, parents=True)

        command = [
            str(a)
            for a in [
                "skewer",
                "-r",
                config.max_error_rate,
                "-d",
                config.max_indel_rate,
                "-m",
                config.mode.value,
                "-l",
                config.min_length,
                "-q",
                config.end_quality,
                "-Q",
                config.mean_quality,
                "-t",
                proc,
                # Skewer spams the console with progress updates. Set quiet to avoid.
                "--quiet",
                # Compress the trimmed output.
                "-z",
                "-o",
                f"{temp_path}/reads",
                *read_paths,
            ]
        ]

        process = await run_subprocess(
            command,
            cwd=read_paths[0].parent,
            env={**os.environ, "LD_LIBRARY_PATH": "/usr/lib/x86_64-linux-gnu"},
        )

        read_paths = await asyncio.to_thread(
            _rename_trimming_results,
            temp_path,
            output_path,
        )

        return SkewerResult(command, output_path, process, read_paths)

    return func


def _rename_trimming_results(temp_path: Path, output_path: Path) -> ReadPaths:
    """Rename Skewer output to a simple name used in Virtool.

    :param path: The path containing the results from Skewer
    """
    shutil.move(
        temp_path / "reads-trimmed.log",
        output_path / "trim.log",
    )

    try:
        return (
            shutil.move(
                temp_path / "reads-trimmed.fastq.gz",
                output_path / "reads_1.fq.gz",
            ),
        )
    except FileNotFoundError:
        return (
            shutil.move(
                temp_path / "reads-trimmed-pair1.fastq.gz",
                output_path / "reads_1.fq.gz",
            ),
            shutil.move(
                temp_path / "reads-trimmed-pair2.fastq.gz",
                output_path / "reads_2.fq.gz",
            ),
        )


def calculate_trimming_min_length(sample: WFSample) -> int:
    """Calculate the minimum trimming length that should be used for the passed sample.

    This takes into account the library type (:class:`.LibraryType`) and the maximum
    observed read length in the sample.

    :param sample: the sample
    :return: the minimum allowed trimmed read length
    """
    if sample.library_type == LibraryType.amplicon:
        return round(0.95 * sample.max_length)

    if sample.max_length < 80:
        return 35

    if sample.max_length < 160:
        return 100

    return 160


async def get_tool_version(
    tool_name: str,
    version_command: list[str],
    run_subprocess: RunSubprocess,
) -> str:
    """Return a parsed tool version from a subprocess version command."""
    output = []

    async def collect_output(line: bytes) -> None:
        output.append(line.decode())

    await run_subprocess(
        version_command,
        stderr_handler=collect_output,
        stdout_handler=collect_output,
    )

    output_text = "".join(output)

    match = re.search(r"\bversion\s+([^\s]+)", output_text)

    if match is None:
        match = re.search(r"\bv?([0-9]+(?:\.[0-9A-Za-z_-]+)+)", output_text)

    if match is None:
        raise ValueError(f"Could not parse {tool_name} version")

    return match.group(1)


async def get_bowtie2_build_version(run_subprocess: RunSubprocess) -> str:
    """Return the version reported by bowtie2-build."""
    return await get_tool_version(
        BOWTIE2_BUILD_TOOL,
        [BOWTIE2_BUILD_TOOL, "--version"],
        run_subprocess,
    )


async def get_skewer_version(run_subprocess: RunSubprocess) -> str:
    """Return the version reported by Skewer."""
    return await get_tool_version(
        SKEWER_TOOL,
        [SKEWER_TOOL, "--version"],
        run_subprocess,
    )


async def get_trimmed_reads_cache_params(
    config: SkewerConfiguration,
    sample: WFSample,
    run_subprocess: RunSubprocess,
) -> dict[str, str | int | float | bool]:
    """Return cache key params for NuVs trimmed reads."""
    return {
        "kind": "trimmed_reads",
        "workflow": WORKFLOW_NAME,
        "workflow_version": WORKFLOW_VERSION,
        "parent_id": sample.id,
        "tool_name": SKEWER_TOOL,
        "tool_version": await get_skewer_version(run_subprocess),
        "min_length": config.min_length,
        "mode": config.mode.value,
        "end_quality": config.end_quality,
        "max_error_rate": config.max_error_rate,
        "max_indel_rate": config.max_indel_rate,
        "mean_quality": config.mean_quality,
    }


async def get_mapping_index_cache_params(
    index_kind: str,
    parent_id: str,
    run_subprocess: RunSubprocess,
    extra_params: dict[str, str] | None = None,
) -> dict[str, str]:
    """Return cache key params for a Bowtie2 mapping index."""
    params = {
        "index_kind": index_kind,
        "workflow": WORKFLOW_NAME,
        "workflow_version": WORKFLOW_VERSION,
        "parent_id": parent_id,
        "tool_name": BOWTIE2_BUILD_TOOL,
        "tool_version": await get_bowtie2_build_version(run_subprocess),
    }

    if extra_params is not None:
        params.update(extra_params)

    return params


def derive_cache_key(params: dict) -> str:
    """Derive a shared workflow cache key from cache params."""
    return derive_key(params)


async def build_bowtie2_index(
    fasta_path: Path,
    index_prefix: Path,
    proc: int,
    run_subprocess: RunSubprocess,
) -> None:
    """Build a Bowtie2 index with a stable local prefix."""
    await run_subprocess(
        [
            BOWTIE2_BUILD_TOOL,
            "--threads",
            str(proc),
            str(fasta_path),
            str(index_prefix),
        ],
    )


async def create_mapping_index(
    cache: WorkflowCache,
    proc: int,
    run_subprocess: RunSubprocess,
    *,
    index_kind: str,
    index_prefix: Path,
    fasta_path: Path,
    parent_id: str,
    extra_params: dict[str, str] | None = None,
) -> None:
    """Build or restore a Bowtie2 mapping index through the shared workflow cache."""
    index_dir = index_prefix.parent
    cache_restore_parent = index_dir.parent
    params = await get_mapping_index_cache_params(
        index_kind,
        parent_id,
        run_subprocess,
        extra_params,
    )
    key = derive_cache_key(params)
    log = logger.bind(
        index_kind=index_kind,
        key=key,
        parent_id=parent_id,
        workflow=WORKFLOW_NAME,
    )

    log.info("checking workflow cache")

    result = await cache.get(key, cache_restore_parent)

    if isinstance(result, CacheHit):
        log.info("restored cached mapping index", outcome="hit")
        return

    log.info("building mapping index", outcome="miss")

    index_dir.mkdir(parents=True, exist_ok=True)

    await build_bowtie2_index(fasta_path, index_prefix, proc, run_subprocess)

    created = await cache.put(key, index_dir, params=params)

    if created:
        log.info("cached built mapping index", outcome="put")
    else:
        log.info("mapping index cache already exists", outcome="put_skipped")


def read_fastq_headers(path: Path) -> set[str]:
    """Return a list of FASTQ headers for the FASTQ file located at `path`.
    Only accepts uncompressed FASTQ files.

    :param path: the path to the FASTQ file
    :return: a list of FASTQ headers

    """
    return {record.id for record in SeqIO.parse(path, "fastq")}


def filter_reads_by_headers(
    headers: set[str],
    out_paths: tuple[Path, Path],
    read_paths: ReadPaths,
):
    """Filter FASTQ data based on whether their headers are in `headers`.

    :param headers: headers for which FASTQ entries should be retained
    :param out_paths: paths to which to write the filtered FASTQ files
    :param read_paths: paths from which to read FASTQ files

    """
    out_handles = tuple(path.open("w") for path in out_paths)

    for path, handle in zip(read_paths, out_handles, strict=False):
        records = (
            record for record in SeqIO.parse(path, "fastq") if record.id in headers
        )

        while True:
            batch = list(itertools.islice(records, 1000))

            if not batch:
                break

            SeqIO.write(batch, handle, "fastq")
