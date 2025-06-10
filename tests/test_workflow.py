import gzip
import shutil
from pathlib import Path

import pytest
from Bio import SeqIO
from syrupy import SnapshotAssertion
from virtool.hmm.models import HMM
from virtool.models.enums import LibraryType
from virtool.samples.models import Sample
from virtool_workflow.analysis.skewer import skewer
from virtool_workflow.analysis.utils import ReadPaths
from virtool_workflow.data.analyses import WFAnalysis
from virtool_workflow.data.hmms import WFHMMs
from virtool_workflow.data.indexes import WFIndex
from virtool_workflow.data.samples import WFSample
from virtool_workflow.data.subtractions import WFSubtraction
from virtool_workflow.pytest_plugin import Data
from virtool_workflow.runtime.run_subprocess import RunSubprocess

from workflow import (
    assemble,
    eliminate_otus,
    eliminate_subtraction,
    process_assembly,
    reunite_pairs,
    trim_reads,
    vfam,
)


@pytest.fixture()
async def analysis(
    data: Data,
    mocker,
) -> WFAnalysis:
    analysis_ = mocker.Mock(spec=WFAnalysis)

    analysis_.id = data.analysis.id
    analysis_.ready = True
    analysis_.files = []
    analysis_.index = data.analysis.index
    analysis_.ml = None
    analysis_.reference = data.analysis.reference
    analysis_.sample = data.analysis.sample
    analysis_.subtractions = data.analysis.subtractions
    analysis_.workflow = "nuvs"

    return analysis_


@pytest.fixture()
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


@pytest.fixture()
def index(data: Data, example_path: Path, work_path: Path) -> WFIndex:
    index_path = work_path / "indexes" / data.index.id
    index_path.parent.mkdir(parents=True)

    shutil.copytree(example_path / "index", index_path)

    return WFIndex(
        id=data.index.id,
        path=index_path,
        manifest={},
        reference=data.index.reference,
        sequence_lengths={},
        sequence_otu_map={},
    )


@pytest.fixture()
async def sample(data: Data, example_path: Path, work_path: Path) -> WFSample:
    path = work_path / "samples" / data.sample.id
    path.mkdir(parents=True)

    shutil.copyfile(example_path / "sample" / "reads_1.fq.gz", path / "reads_1.fq.gz")

    return WFSample(
        id=data.sample.id,
        library_type=LibraryType.normal,
        name=data.sample.name,
        paired=False,
        quality=data.sample.quality,
        read_paths=(path / "reads_1.fq.gz",),
    )


@pytest.fixture()
async def subtractions(
    data: Data,
    example_path: Path,
    work_path: Path,
) -> list[WFSubtraction]:
    subtractions_path = work_path / "subtractions"
    subtractions_path.mkdir(parents=True)

    subtraction_1, subtraction_2 = (
        WFSubtraction(
            id=f"subtraction_{suffix}",
            files=[],
            gc=data.subtraction.gc,
            nickname=f"Subby {suffix}",
            name=f"Subtraction {suffix}",
            path=subtractions_path / f"subtraction_{suffix}",
        )
        for suffix in (1, 2)
    )

    for subtraction in (subtraction_1, subtraction_2):
        shutil.copytree(example_path / "subtraction", subtraction.path)

    return [subtraction_1, subtraction_2]


@pytest.fixture()
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
    work_path: Path,
):
    proc = 2

    sample.quality.length = [45, 76]

    # Instantiate the skewer fixture.
    skewer_ = skewer(2, run_subprocess)

    await trim_reads(proc, sample, skewer_, work_path)

    with gzip.open(work_path / "trimmed" / "reads_1.fq.gz", "rt") as f:
        assert {
            (record.id, record.seq)
            for record in SeqIO.parse(
                f,
                "fastq",
            )
        } == snapshot


async def test_eliminate_otus(
    index: WFIndex,
    run_subprocess: RunSubprocess,
    sample: WFSample,
    snapshot: SnapshotAssertion,
    trimmed_read_paths: ReadPaths,
    work_path: Path,
):
    """Make sure the eliminated output file, unmapped_otus.fq, matches the expected."""
    proc = 1

    await eliminate_otus(index, proc, run_subprocess, trimmed_read_paths, work_path)

    with open(work_path / "unmapped_otus.fq") as f:
        assert [line.rstrip() for line in f] == snapshot


@pytest.mark.parametrize("no_subtractions", [True, False])
async def test_eliminate_subtraction(
    no_subtractions: bool,
    example_path: Path,
    run_subprocess: RunSubprocess,
    snapshot: SnapshotAssertion,
    subtractions: list[WFSubtraction],
    work_path,
):
    """Test that the step eliminates subractions when provided, but does nothing if no
    subtractions are provided.
    """
    if no_subtractions:
        subtractions = []

    shutil.copyfile(example_path / "unmapped_otus.fq", work_path / "unmapped_otus.fq")

    await eliminate_subtraction(2, run_subprocess, subtractions, work_path)

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
    sample: Sample,
    snapshot: SnapshotAssertion,
    run_subprocess: RunSubprocess,
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
