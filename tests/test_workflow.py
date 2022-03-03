import gzip
import json
import shutil
from pathlib import Path
from shutil import copytree
from typing import List

import pytest
import virtool_workflow.execution.run_subprocess
from Bio import SeqIO
from aiohttp.test_utils import make_mocked_coro
from virtool_workflow.analysis.analysis import Analysis
from virtool_workflow.analysis.hmms import HMMs
from virtool_workflow.analysis.indexes import Index
from virtool_workflow.analysis.library_types import LibraryType
from virtool_workflow.analysis.reads import Reads
from virtool_workflow.data_model import (
    NucleotideComposition,
    Subtraction,
    Reference,
    Sample,
    HMM,
)

from workflow import (
    eliminate_otus,
    eliminate_subtraction,
    reunite_pairs,
    assemble,
    process_fasta,
    vfam,
)

TEST_DATA_PATH = Path(__file__).parent / "data"
FASTQ_PATH = TEST_DATA_PATH / "test.fq"
INDEX_PATH = TEST_DATA_PATH / "index"
SUBTRACTION_PATH = TEST_DATA_PATH / "subtraction"


@pytest.fixture
def work_path(tmpdir):
    return Path(tmpdir) / "work"


@pytest.fixture
def run_in_executor():
    async def _run_in_executor(func, *args):
        return func(*args)

    return _run_in_executor


@pytest.fixture
def run_subprocess():
    return virtool_workflow.execution.run_subprocess.run_subprocess()


@pytest.fixture
def hmms(work_path: Path):
    hmms_path = work_path / "hmms"

    shutil.copytree(TEST_DATA_PATH / "hmms", hmms_path)

    annotations = [
        HMM(
            id="foo",
            cluster=2,
            count=21,
            entries=list(),
            families=dict(),
            genera=dict(),
            hidden=False,
            length=42,
            mean_entropy=0.0001,
            total_entropy=0.0001,
            names=("Test", "Foo", "Bar"),
        ),
        HMM(
            id="bar",
            cluster=9,
            count=21,
            entries=list(),
            families=dict(),
            genera=dict(),
            hidden=False,
            length=42,
            mean_entropy=0.0001,
            total_entropy=0.0001,
            names=("Test", "Foo", "Bar"),
        ),
    ]

    return HMMs(annotations, hmms_path)


@pytest.fixture
def indexes(run_in_executor, run_subprocess, work_path) -> List[Index]:
    index_path = work_path / "references"
    index_path.mkdir(parents=True)

    index_path = index_path / "foo"

    shutil.copytree(INDEX_PATH, index_path)

    reference = Reference(
        id="reference_1",
        data_type="genome",
        description="Reference 1",
        name="Reference 1",
        organism="viruses",
    )

    index = Index(
        "foo", dict(), reference, True, index_path, run_in_executor, run_subprocess
    )

    return [index]


@pytest.fixture
async def sample():
    return Sample(
        id="sample",
        name="Sample 1",
        host="",
        isolate="",
        locale="",
        library_type=LibraryType.other,
        paired=False,
        quality=dict(),
        nuvs=False,
    )


@pytest.fixture
async def reads(sample: Sample, work_path: Path):
    reads_path = work_path / "reads"
    reads_path.mkdir(parents=True)

    shutil.copy(FASTQ_PATH, reads_path / "reads_1.fq.gz")

    return Reads(sample, {}, reads_path)


@pytest.fixture
async def analysis(
    indexes: List[Index], sample: Sample, subtractions: List[Subtraction]
):
    upload_files = make_mocked_coro()

    return Analysis(
        upload_files,
        id="foo",
        files=[],
        sample=sample,
        index=indexes[0],
        subtractions=subtractions,
        ready=True,
    )


@pytest.fixture
async def subtractions(work_path):
    subtractions_path = work_path / "subtractions"
    subtractions_path.mkdir(parents=True)

    subtraction_path = work_path / "subtractions" / "subtraction"

    copytree(SUBTRACTION_PATH, subtraction_path)

    nucleotide_composition = NucleotideComposition(
        a=0.1,
        t=0.2,
        g=0.3,
        c=0.4,
    )

    subtraction = Subtraction(
        id="arabidopsis_thaliana",
        name="Arabidopsis thaliana",
        nickname="Thalecress",
        count=12,
        gc=nucleotide_composition,
        path=subtraction_path,
    )

    return [subtraction]


async def test_eliminate_otus(indexes, reads: Reads, run_subprocess, work_path: Path):
    await eliminate_otus(indexes, 2, reads, run_subprocess, work_path)

    actual_path = work_path / "unmapped_otus.fq"

    with open(actual_path, "r") as f:
        actual = [line.rstrip() for line in f]
        actual = {tuple(actual[i : i + 4]) for i in range(0, len(actual), 4)}

    assert len(actual) > 0

    with open(TEST_DATA_PATH / "unmapped_otus.fq", "r") as f:
        expected = [line.rstrip() for line in f]
        expected = {tuple(expected[i : i + 4]) for i in range(0, len(expected), 4)}

    assert actual == expected


async def test_eliminate_subtraction(run_subprocess, subtractions, work_path):
    shutil.copy(TEST_DATA_PATH / "unmapped_otus.fq", work_path / "unmapped_otus.fq")
    await eliminate_subtraction(2, run_subprocess, subtractions, work_path)


@pytest.mark.parametrize("paired", [False, True], ids=["unpaired", "paired"])
async def test_reunite_pairs(paired, run_in_executor, reads, sample, work_path):
    if paired:
        sample.paired = paired

        unite_path = TEST_DATA_PATH / "unite.json"

        with open(unite_path, "r") as f:
            unite = json.load(f)

        for path, key in [(reads.left, "left"), (reads.right, "right")]:
            with gzip.open(path, "wt") as f:
                for line in unite[key]:
                    f.write(line + "\n")

        separate_path = work_path / "unmapped_hosts.fq"

        with open(separate_path, "w") as f:
            for line in unite["separate"]:
                f.write(line + "\n")

    await reunite_pairs(run_in_executor, reads, work_path)

    if paired:
        for filename, key in [
            ("unmapped_1.fq", "united_left"),
            ("unmapped_2.fq", "united_right"),
        ]:
            with open(work_path / filename, "r") as f:
                assert [line.rstrip() for line in f] == unite[key]


@pytest.mark.parametrize("paired", [False, True], ids=["unpaired", "paired"])
async def test_assemble(
    paired: bool,
    analysis: Analysis,
    sample: Sample,
    run_in_executor,
    run_subprocess,
    work_path: Path,
):
    sample.paired = paired

    proc = 2
    mem = 10

    if paired:
        for suffix in (1, 2):
            shutil.copy(
                TEST_DATA_PATH / "unmapped_{}.fq".format(suffix),
                work_path / "unmapped_{}.fq".format(suffix),
            )
    else:
        shutil.copy(TEST_DATA_PATH / "unmapped_1.fq", work_path / "unmapped_hosts.fq")

    await assemble(
        analysis, mem, proc, run_in_executor, run_subprocess, sample, work_path
    )

    expected_path = TEST_DATA_PATH / "scaffolds_{}.fa".format("p" if paired else "u")
    assembly_path = work_path / "assembly.fa"
    compressed_path = work_path / "assembly.fa.gz"

    expected = {
        (record.id, record.seq) for record in SeqIO.parse(expected_path, "fasta")
    }

    # Check that scaffolds.fasta from SPAdes matches expected data.
    assert {
        (record.id, record.seq) for record in SeqIO.parse(assembly_path, "fasta")
    } == expected

    # Check that compressed assembly matches expected data.
    with gzip.open(compressed_path, "rt") as f:
        assert {
            (record.id, record.seq) for record in SeqIO.parse(f, "fasta")
        } == expected

    assert analysis.to_upload == [(compressed_path, "fasta")]


async def test_process_fasta(
    data_regression,
    file_regression,
    analysis: Analysis,
    run_in_executor,
    work_path: Path,
):
    spades_path = work_path / "spades"
    spades_path.mkdir(parents=True)

    shutil.copy(TEST_DATA_PATH / "scaffolds_u.fa", work_path / "assembly.fa")

    results = dict()

    await process_fasta(analysis, 2, results, run_in_executor, work_path)

    data_regression.check(results)

    with open(work_path / "orfs.fa") as f:
        file_regression.check(f.read())


async def test_vfam(
    data_regression, analysis: Analysis, hmms: HMMs, run_subprocess, work_path: Path
):
    with open(work_path / "orfs.fa", "w") as f:
        f.write(">sequence_0.0\n")
        f.write(
            "MVAVRAPRRKRASATDLYKTCKAAGTCPPDVIPKIEGSTLADKILQWSGLGIFLGGLGIGTGTGSGGRTGYIPLGGGGRPSVVDIGPTRPPIIIEPVGPTEPSIVT"
            "LVEESSIIQSGAPIPTFSGGNGFELTTSSATTPAVLDITPSAGTVHVTSTNIQNPLYIEPPIDIPQAGEASGHIFTTTSTAGTHSYEEIPMEVFASTNGTGLEPIS"
            "STPIPGIQRVSAPRLYSKAYQQVKVTDPNFIGNPSTFVTFDNPAYEPIDETLTYASSSTVAPDPDFLDIIALHRPALTSRKGTVRYSRLGQKATMKTRSGKQIGAT"
            "VHYYHDISPIQSFAEHEEIELQPLHTSTHSSAPLFDIYADPDTVPSIHTPRMSYSPTTLPVPRYASNVFSSINTSTTNVTVPLSTSFELPVYSGSDIYTPTSSPTW"
            "PSLPPPPTTNLPAIVVHGDNYYLWPYIYLIHKRRKRMPYFFSDGFVAY"
        )
        f.write("\n>sequence_2.0\n")
        f.write(
            "MSSLVSETSNSEVGSQMESPGRGGQSIDAPSSSCFKVRARNLFLTYSKCNLTAVFLLEYISSLLKKYCPTYIYVAQEAHKDGSHHLHCIIQCSKYVRTTSAKFFDI"
            "KEFHPNVQNPRMPKKALSYCKKSPISEAEYGVFQEIKRPRKKKADAPSTKDAKMAEIIKSSTNKEDYLSMVRKSFPFDWATRLQQFQFSAESLFPSTPPPYVDPFG"
            "MPSQDTHPVIGAWLRDELYTDRSPTERRRSLYICGPTRTGKTSWARSLGSHNYWQHSVDFLHVIQNARYNVIDDIPFKFVPCWKGLVGSQKDITVNPKYGKKRLLS"
            "NGIPCIILVNEDEDWLQQMQPSQADWFNANAVVHYMYSGESFFEAL"
        )

    results = {
        "hits": [
            {"orfs": [{"name": "Foo"}]},
            {"orfs": [{"name": "Nil"}]},
            {"orfs": [{"name": "Bar"}]},
        ]
    }

    await vfam(analysis, hmms, 2, results, run_subprocess, work_path)

    data_regression.check(results)
