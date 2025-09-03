import asyncio
import itertools
import os
import shutil
from asyncio.subprocess import Process
from dataclasses import dataclass
from enum import Enum
from pathlib import Path
from tempfile import mkdtemp
from typing import Protocol

from Bio import SeqIO
from pyfixtures import fixture
from virtool.models.enums import LibraryType
from virtool.workflow import RunSubprocess
from virtool.workflow.analysis import ReadPaths
from virtool.workflow.data.samples import WFSample


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
    """Represents the result of running Skewer to trim a paired or unpaired FASTQ dataset."""

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
        """The path to one of:

        - the FASTQ trimming result for an unpaired Illumina dataset
        - the FASTA trimming result for the left reads of a paired Illumina dataset

        """
        return self.read_paths[0]

    @property
    def right(self) -> Path | None:
        """The path to the rights reads of a paired Illumina dataset.

        ``None`` if the dataset in unpaired.

        :type: :class:`.Path`

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


def calculate_trimming_cache_key(
    sample_id: str,
    trimming_parameters: dict,
    program: str = "skewer",
):
    """Compute a unique cache key.

    **This is not currently used.**

    :param sample_id: The ID of the sample being trimmed.
    :param trimming_parameters: The trimming parameters.
    :param program: The name of the trimming program.
    :return: A unique cache key.

    """
    raw_key = "reads-" + json.dumps(
        {
            "id": sample_id,
            "parameters": trimming_parameters,
            "program": program,
        },
        sort_keys=True,
    )

    return hashlib.sha256(raw_key.encode()).hexdigest()


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
    """Filter paired, input FASTQ files based on whether their headers are in
    `headers`.
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
