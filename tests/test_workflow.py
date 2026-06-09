import asyncio
import gzip
import shutil
from collections.abc import Callable
from pathlib import Path

import pytest
from Bio import SeqIO
from pytest_mock import MockerFixture
from structlog import get_logger
from syrupy import SnapshotAssertion
from virtool.caches.utils import derive_key
from virtool.hmm.models import HMM
from virtool.models.enums import LibraryType
from virtool.samples.models import Sample
from virtool.workflow.analysis import ReadPaths
from virtool.workflow.data.analyses import WFAnalysis
from virtool.workflow.data.cache import WorkflowCache
from virtool.workflow.data.hmms import WFHMMs
from virtool.workflow.data.indexes import WFIndex
from virtool.workflow.data.samples import WFSample
from virtool.workflow.data.subtractions import WFSubtraction
from virtool.workflow.errors import JobsAPINotFoundError
from virtool.workflow.pytest_plugin import WorkflowData
from virtool.workflow.runtime.run_subprocess import RunSubprocess

from utils import (
    SkewerConfiguration,
    SkewerMode,
    get_mapping_index_cache_params,
    get_trimmed_reads_cache_params,
    skewer,
)
from workflow import (
    assemble,
    create_reference_index,
    create_subtraction_indexes,
    eliminate_otus,
    eliminate_subtraction,
    get_subtraction_index_path,
    process_assembly,
    reunite_pairs,
    trim_reads,
    vfam,
)

logger = get_logger("test")

BOWTIE2_INDEX_SUFFIXES = (
    "1.bt2",
    "2.bt2",
    "3.bt2",
    "4.bt2",
    "rev.1.bt2",
    "rev.2.bt2",
)

SubtractionIndexPath = Callable[[int | WFSubtraction], Path]


class _FakeWorkflowCacheAPI:
    """Fake only the API calls used by the real workflow cache."""

    def __init__(self, work_dir: Path) -> None:
        self.work_dir = work_dir
        self.stored: dict[str, tuple[Path, dict | None]] = {}

    async def get_cache(self, key: str, dest: Path) -> None:
        try:
            source, _ = self.stored[key]
        except KeyError:
            raise JobsAPINotFoundError from None

        dest.parent.mkdir(parents=True, exist_ok=True)
        await asyncio.to_thread(shutil.copyfile, source, dest)

    async def put_cache(
        self,
        key: str,
        path: Path,
        params: dict | None = None,
    ) -> bool:
        if key in self.stored:
            return False

        stored_path = self.work_dir / "stored" / key / "cache.tar"
        stored_path.parent.mkdir(parents=True, exist_ok=True)
        await asyncio.to_thread(shutil.copyfile, path, stored_path)
        self.stored[key] = (stored_path, params)

        return True


@pytest.fixture
def workflow_cache(tmp_path: Path) -> WorkflowCache:
    api = _FakeWorkflowCacheAPI(tmp_path / "fake_workflow_cache_api")
    return WorkflowCache(api, tmp_path / "workflow_cache")


def write_bowtie2_bundle(path: Path, prefix: str, content: bytes = b"cached"):
    path.mkdir(parents=True, exist_ok=True)

    for suffix in BOWTIE2_INDEX_SUFFIXES:
        (path / f"{prefix}.{suffix}").write_bytes(content)


def read_directory_bytes(path: Path) -> dict[str, bytes]:
    return {
        child.name: child.read_bytes() for child in path.iterdir() if child.is_file()
    }


def bowtie2_bundle_bytes(prefix: str, content: bytes | None = None) -> dict[str, bytes]:
    return {
        f"{prefix}.{suffix}": content or f"{prefix}.{suffix}".encode()
        for suffix in BOWTIE2_INDEX_SUFFIXES
    }


def assert_bowtie2_bundle_exists(prefix: Path) -> None:
    """Assert that Bowtie2 wrote the expected index files for a prefix."""
    for suffix in BOWTIE2_INDEX_SUFFIXES:
        path = prefix.parent / f"{prefix.name}.{suffix}"
        assert path.exists()
        assert path.stat().st_size > 0


@pytest.fixture
async def analysis(
    mocker: MockerFixture,
    workflow_data: WorkflowData,
) -> WFAnalysis:
    analysis_ = mocker.Mock(spec=WFAnalysis)

    analysis_.id = workflow_data.analysis.id
    analysis_.ready = True
    analysis_.files = []
    analysis_.index = workflow_data.analysis.index
    analysis_.ml = None
    analysis_.reference = workflow_data.analysis.reference
    analysis_.sample = workflow_data.analysis.sample
    analysis_.subtractions = workflow_data.analysis.subtractions
    analysis_.workflow = "nuvs"

    return analysis_


@pytest.fixture
def hmms(example_path: Path, work_path: Path) -> WFHMMs:
    hmms_path = work_path / "hmms"

    shutil.copytree(example_path / "hmms", hmms_path)

    annotations = [
        HMM(
            id="foo",
            cluster=2,
            count=21,
            entries=[],
            families={},
            genera={},
            length=42,
            mean_entropy=0.0001,
            total_entropy=0.0001,
            names=["Test", "Foo", "Bar"],
        ),
        HMM(
            id="bar",
            cluster=9,
            count=21,
            entries=[],
            families={},
            genera={},
            length=42,
            mean_entropy=0.0001,
            total_entropy=0.0001,
            names=["Test", "Foo", "Bar"],
        ),
    ]

    return WFHMMs(annotations, hmms_path)


@pytest.fixture
def index(workflow_data: WorkflowData, example_path: Path, work_path: Path) -> WFIndex:
    index_path = work_path / "indexes" / workflow_data.index.id
    index_path.parent.mkdir(parents=True)

    shutil.copytree(example_path / "index", index_path)
    index = WFIndex(
        id=workflow_data.index.id,
        path=index_path,
        manifest={},
        reference=workflow_data.index.reference,
        sequence_lengths={},
        sequence_otu_map={},
    )

    index.fasta_path.write_text(">reference\nACGT\n")

    return index


@pytest.fixture
def reference_index_path(work_path: Path) -> Path:
    return work_path / "reference_index" / "reference"


@pytest.fixture
async def sample(
    workflow_data: WorkflowData,
    example_path: Path,
    work_path: Path,
) -> WFSample:
    path = work_path / "samples" / workflow_data.sample.id
    path.mkdir(parents=True)

    shutil.copyfile(example_path / "sample" / "reads_1.fq.gz", path / "reads_1.fq.gz")

    return WFSample(
        id=workflow_data.sample.id,
        library_type=LibraryType.normal,
        name=workflow_data.sample.name,
        paired=False,
        quality=workflow_data.sample.quality,
        read_paths=(path / "reads_1.fq.gz",),
    )


@pytest.fixture
async def subtractions(
    example_path: Path,
    workflow_data: WorkflowData,
    work_path: Path,
) -> list[WFSubtraction]:
    subtractions_path = work_path / "subtractions"
    subtractions_path.mkdir(parents=True)

    subtraction_1, subtraction_2 = (
        WFSubtraction(
            id=f"subtraction_{suffix}",
            files=[],
            gc=workflow_data.subtraction.gc,
            nickname=f"Subby {suffix}",
            name=f"Subtraction {suffix}",
            path=subtractions_path / f"subtraction_{suffix}",
        )
        for suffix in (1, 2)
    )

    for subtraction in (subtraction_1, subtraction_2):
        shutil.copytree(
            example_path / "subtractions" / "arabidopsis_thaliana",
            subtraction.path,
        )

    return [subtraction_1, subtraction_2]


@pytest.fixture
def subtraction_indexes_path(work_path: Path) -> Path:
    return work_path / "subtraction_indexes"


@pytest.fixture
def subtraction_index_path(
    subtractions: list[WFSubtraction],
    subtraction_indexes_path: Path,
) -> SubtractionIndexPath:
    def _subtraction_index_path(subtraction: int | WFSubtraction = 0) -> Path:
        if isinstance(subtraction, int):
            subtraction = subtractions[subtraction]

        return get_subtraction_index_path(subtraction_indexes_path, subtraction.id)

    return _subtraction_index_path


@pytest.fixture
def trimmed_read_paths(example_path: Path, work_path: Path) -> ReadPaths:
    path = work_path / "trimmed"
    path.mkdir()

    left_path = path / "reads_1.fq.gz"

    shutil.copyfile(example_path / "sample" / "reads_1.fq.gz", left_path)

    return (left_path,)


async def test_trim_reads(
    run_subprocess: RunSubprocess,
    sample: WFSample,
    snapshot: SnapshotAssertion,
    workflow_cache: WorkflowCache,
    work_path: Path,
):
    proc = 2

    sample.quality.length = [45, 76]

    # Instantiate the skewer fixture.
    skewer_ = skewer(2, run_subprocess)
    await trim_reads(
        workflow_cache,
        proc,
        run_subprocess,
        sample,
        skewer_,
        work_path,
    )

    with gzip.open(work_path / "trimmed" / "reads_1.fq.gz", "rt") as f:
        assert {
            (record.id, record.seq)
            for record in SeqIO.parse(
                f,
                "fastq",
            )
        } == snapshot


async def test_get_trimmed_reads_cache_params(
    run_subprocess: RunSubprocess,
    sample: WFSample,
):
    sample.paired = True
    sample.library_type = LibraryType.normal
    config = SkewerConfiguration(
        min_length=35,
        mode=SkewerMode.PAIRED_END,
        number_of_processes=4,
    )

    params = await get_trimmed_reads_cache_params(config, sample, run_subprocess)

    assert (
        params.items()
        >= {
            "kind": "trimmed_reads",
            "workflow": "nuvs",
            "parent_id": sample.id,
            "tool_name": "skewer",
            "tool_version": "0.2.2",
            "min_length": 35,
            "mode": "pe",
            "end_quality": 20,
            "max_error_rate": 0.1,
            "max_indel_rate": 0.03,
            "mean_quality": 25,
        }.items()
    )


@pytest.mark.parametrize("paired", [False, True], ids=["unpaired", "paired"])
async def test_trim_reads_cache_hit(
    paired: bool,
    run_subprocess: RunSubprocess,
    sample: WFSample,
    tmp_path: Path,
    workflow_cache: WorkflowCache,
    work_path: Path,
):
    sample.paired = paired
    sample.quality.length = [45, 76]

    if paired:
        right_path = tmp_path / "reads_2.fq.gz"
        shutil.copyfile(sample.read_paths[0], right_path)
        sample.read_paths = (*sample.read_paths, right_path)

    cached_path = tmp_path / "trimmed"
    cached_path.mkdir()
    (cached_path / "reads_1.fq.gz").write_bytes(b"cached-1")

    if paired:
        (cached_path / "reads_2.fq.gz").write_bytes(b"cached-2")

    skewer_ = skewer(4, run_subprocess)

    config = SkewerConfiguration(
        min_length=35,
        mode=SkewerMode.PAIRED_END if paired else SkewerMode.SINGLE_END,
        number_of_processes=4,
    )
    params = await get_trimmed_reads_cache_params(config, sample, run_subprocess)
    key = derive_key(params)

    assert await workflow_cache.put(key, cached_path, params)

    await trim_reads(workflow_cache, 4, run_subprocess, sample, skewer_, work_path)

    assert (work_path / "trimmed" / "reads_1.fq.gz").read_bytes() == b"cached-1"

    if paired:
        assert (work_path / "trimmed" / "reads_2.fq.gz").read_bytes() == b"cached-2"


@pytest.mark.parametrize("paired", [False, True], ids=["unpaired", "paired"])
async def test_trim_reads_cache_miss(
    paired: bool,
    run_subprocess: RunSubprocess,
    sample: WFSample,
    tmp_path: Path,
    workflow_cache: WorkflowCache,
    work_path: Path,
):
    sample.paired = paired
    sample.quality.length = [45, 76]

    if paired:
        right_path = tmp_path / "reads_2.fq.gz"
        shutil.copyfile(sample.read_paths[0], right_path)
        sample.read_paths = (*sample.read_paths, right_path)

    skewer_ = skewer(4, run_subprocess)

    await trim_reads(workflow_cache, 4, run_subprocess, sample, skewer_, work_path)

    assert (work_path / "trimmed" / "reads_1.fq.gz").stat().st_size > 0

    if paired:
        assert (work_path / "trimmed" / "reads_2.fq.gz").stat().st_size > 0


async def test_create_reference_index_cache_hit(
    index: WFIndex,
    reference_index_path: Path,
    run_subprocess: RunSubprocess,
    tmp_path: Path,
    workflow_cache: WorkflowCache,
):
    source = tmp_path / reference_index_path.parent.name
    write_bowtie2_bundle(source, "reference")

    params = await get_mapping_index_cache_params(
        "reference_mapping_index",
        index.id,
        run_subprocess,
        {"source": "reference_fasta"},
    )
    key = derive_key(params)

    assert await workflow_cache.put(key, source, params)

    result = await create_reference_index(
        workflow_cache,
        index,
        4,
        reference_index_path,
        run_subprocess,
    )

    assert result == reference_index_path
    assert read_directory_bytes(reference_index_path.parent) == read_directory_bytes(
        source,
    )


async def test_create_reference_index_cache_miss(
    index: WFIndex,
    reference_index_path: Path,
    run_subprocess: RunSubprocess,
    workflow_cache: WorkflowCache,
):
    result = await create_reference_index(
        workflow_cache,
        index,
        4,
        reference_index_path,
        run_subprocess,
    )

    assert result == reference_index_path
    assert_bowtie2_bundle_exists(reference_index_path)


async def test_workflow_cache_put_is_skipped_when_key_is_already_stored(
    index: WFIndex,
    reference_index_path: Path,
    run_subprocess: RunSubprocess,
    workflow_cache: WorkflowCache,
):
    params = await get_mapping_index_cache_params(
        "reference_mapping_index",
        index.id,
        run_subprocess,
        {"source": "reference_fasta"},
    )
    key = derive_key(params)
    path = reference_index_path.parent

    write_bowtie2_bundle(path, "reference")
    assert await workflow_cache.put(key, path, params)
    assert not await workflow_cache.put(key, path, params)


async def test_create_subtraction_indexes_cache_hit(
    run_subprocess: RunSubprocess,
    subtraction_index_path: SubtractionIndexPath,
    subtraction_indexes_path: Path,
    subtractions: list[WFSubtraction],
    tmp_path: Path,
    workflow_cache: WorkflowCache,
):
    subtraction = subtractions[0]
    index_path = subtraction_index_path(0)
    source = tmp_path / index_path.parent.name
    write_bowtie2_bundle(source, "subtraction")

    params = await get_mapping_index_cache_params(
        "subtraction_mapping_index",
        subtraction.id,
        run_subprocess,
        {"source": "subtraction_fasta"},
    )
    key = derive_key(params)

    assert await workflow_cache.put(key, source, params)

    result = await create_subtraction_indexes(
        workflow_cache,
        4,
        run_subprocess,
        subtraction_indexes_path,
        [subtraction],
    )

    assert result is None
    assert read_directory_bytes(index_path.parent) == read_directory_bytes(
        source,
    )


async def test_create_subtraction_indexes_cache_miss(
    run_subprocess: RunSubprocess,
    subtraction_index_path: SubtractionIndexPath,
    subtraction_indexes_path: Path,
    subtractions: list[WFSubtraction],
    workflow_cache: WorkflowCache,
):
    subtraction = subtractions[0]
    index_path = subtraction_index_path(0)

    result = await create_subtraction_indexes(
        workflow_cache,
        4,
        run_subprocess,
        subtraction_indexes_path,
        [subtraction],
    )

    assert result is None
    assert_bowtie2_bundle_exists(index_path)


async def test_eliminate_otus(
    example_path: Path,
    reference_index_path: Path,
    run_subprocess: RunSubprocess,
    snapshot: SnapshotAssertion,
    trimmed_read_paths: ReadPaths,
    work_path: Path,
):
    """Make sure the eliminated output file, unmapped_otus.fq, matches the expected."""
    proc = 1

    source = work_path / reference_index_path.parent.name
    shutil.copytree(example_path / "index", source)

    await eliminate_otus(
        proc,
        reference_index_path,
        run_subprocess,
        trimmed_read_paths,
        work_path,
    )

    with open(work_path / "unmapped_otus.fq") as f:
        assert [line.rstrip() for line in f] == snapshot


@pytest.mark.parametrize("no_subtractions", [True, False])
async def test_eliminate_subtraction(
    no_subtractions: bool,
    example_path: Path,
    run_subprocess: RunSubprocess,
    snapshot: SnapshotAssertion,
    subtraction_index_path: SubtractionIndexPath,
    subtraction_indexes_path: Path,
    subtractions: list[WFSubtraction],
    work_path: Path,
):
    """Test that the step eliminates subractions when provided, but does nothing if no
    subtractions are provided.
    """
    if no_subtractions:
        subtractions = []

    shutil.copyfile(example_path / "unmapped_otus.fq", work_path / "unmapped_otus.fq")

    if subtractions:
        for index in range(len(subtractions)):
            index_path = subtraction_index_path(index)
            shutil.copytree(
                example_path / "subtractions" / "arabidopsis_thaliana",
                index_path.parent,
            )

    await eliminate_subtraction(
        2,
        run_subprocess,
        subtractions,
        subtraction_indexes_path,
        work_path,
    )

    if no_subtractions:
        assert (
            open(work_path / "unmapped_subtractions.fq").read()
            == open(
                work_path / "unmapped_otus.fq",
            ).read()
        )

    assert {
        (record.id, record.seq)
        for record in SeqIO.parse(
            work_path / work_path / "unmapped_subtractions.fq",
            "fastq",
        )
    } == snapshot


@pytest.mark.parametrize("paired", [False, True], ids=["unpaired", "paired"])
async def test_reunite_pairs(
    example_path: Path,
    paired: bool,
    sample: WFSample,
    snapshot: SnapshotAssertion,
    unite: dict[str, list[str]],
    work_path: Path,
):
    trimmed_path = work_path / "trimmed"
    trimmed_path.mkdir()

    trimmed_read_paths = (
        work_path / "trimmed" / "reads_1.fq.gz",
        work_path / "trimmed" / "reads_2.fq.gz",
    )

    if paired:
        sample.paired = paired

        for suffix in (1, 2):
            shutil.copy(
                example_path / "sample" / f"paired_{suffix}.fq.gz",
                work_path / "trimmed" / f"reads_{suffix}.fq.gz",
            )

        for path, key in zip(trimmed_read_paths, ("left", "right"), strict=False):
            with gzip.open(path, "wt") as f:
                for line in unite[key]:
                    f.write(line + "\n")

        with open(work_path / "unmapped_subtractions.fq", "w") as f:
            for line in unite["separate"]:
                f.write(line + "\n")

    await reunite_pairs(
        2,
        sample,
        trimmed_read_paths,
        work_path,
    )

    if paired:
        for filename in ("unmapped_1.fq", "unmapped_2.fq"):
            with open(work_path / filename) as f:
                assert sorted(line.rstrip() for line in f) == snapshot(name=filename)


@pytest.mark.parametrize("paired", [False, True], ids=["unpaired", "paired"])
async def test_assemble(
    paired: bool,
    analysis: WFAnalysis,
    example_path: Path,
    run_subprocess: RunSubprocess,
    sample: Sample,
    snapshot: SnapshotAssertion,
    work_path: Path,
):
    sample.paired = paired

    if paired:
        for suffix in (1, 2):
            filename = f"unmapped_{suffix}.fq"

            shutil.copy(
                example_path / filename,
                work_path / filename,
            )
    else:
        shutil.copy(
            example_path / "unmapped_1.fq",
            work_path / "unmapped_subtractions.fq",
        )

    mem = 5
    proc = 1

    await assemble(analysis, mem, proc, run_subprocess, sample, work_path)

    fasta = {
        (record.id, record.seq)
        for record in SeqIO.parse(work_path / "spades/scaffolds.fasta", "fasta")
    }

    compressed_path = work_path / "assembly.fa.gz"

    # Check that scaffolds.fasta and assembly.fa.gz match expected data. The latter is
    # uploaded to the API.
    with gzip.open(compressed_path, "rt") as f:
        assert (
            {(record.id, record.seq) for record in SeqIO.parse(f, "fasta")}
            == fasta
            == snapshot(name="scaffolds.fasta")
        )

    analysis.upload_file.assert_called_with(compressed_path, "fasta")


async def test_process_fasta(
    analysis: WFAnalysis,
    example_path: Path,
    snapshot: SnapshotAssertion,
    work_path: Path,
):
    spades_path = work_path / "spades"
    spades_path.mkdir(parents=True)

    shutil.copy(example_path / "scaffolds.fasta", spades_path / "scaffolds.fasta")

    results = {}

    await process_assembly(analysis, 2, results, work_path)

    assert results == snapshot(name="results")
    assert open(work_path / "orfs.fa").read() == snapshot(name="fasta")


async def test_vfam(
    analysis: WFAnalysis,
    example_path: Path,
    hmms: WFHMMs,
    run_subprocess: RunSubprocess,
    snapshot: SnapshotAssertion,
    work_path: Path,
):
    """Test that VFAM results match expected and are correctly written to the results
    dictionary.
    """
    shutil.copyfile(example_path / "orfs.fa", work_path / "orfs.fa")

    results = {
        "hits": [
            {"orfs": [{"name": "Foo"}]},
            {"orfs": [{"name": "Nil"}]},
            {"orfs": [{"name": "Bar"}]},
        ],
    }

    await vfam(analysis, hmms, 2, results, run_subprocess, work_path)

    assert results == snapshot(name="results")
