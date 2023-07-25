import gzip
import json
import logging
import shutil
from pathlib import Path
from shutil import copytree
from typing import List

import arrow
import pytest
import virtool_workflow.runtime.run_subprocess
from Bio import SeqIO
from aiohttp.test_utils import make_mocked_coro
from faker import Faker
from pydantic_factories import ModelFactory, Use
from virtool_core.models.enums import LibraryType
from virtool_core.models.hmm import HMM
from virtool_core.models.index import Index
from virtool_core.models.job import JobNested
from virtool_core.models.reference import ReferenceNested, ReferenceDataType
from virtool_core.models.samples import Sample
from virtool_core.models.subtraction import NucleotideComposition, SubtractionUpload
from virtool_core.models.user import UserNested
from virtool_workflow.analysis.hmms import HMMs
from virtool_workflow.analysis.reads import Reads
from virtool_workflow.data_model.analysis import WFAnalysis
from virtool_workflow.data_model.indexes import WFIndex
from virtool_workflow.data_model.samples import WFSample
from virtool_workflow.data_model.subtractions import WFSubtraction

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
def run_subprocess():
    return virtool_workflow.runtime.run_subprocess.run_subprocess()


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
            length=42,
            mean_entropy=0.0001,
            total_entropy=0.0001,
            names=["Test", "Foo", "Bar"],
        ),
        HMM(
            id="bar",
            cluster=9,
            count=21,
            entries=list(),
            families=dict(),
            genera=dict(),
            length=42,
            mean_entropy=0.0001,
            total_entropy=0.0001,
            names=["Test", "Foo", "Bar"],
        ),
    ]

    return HMMs(annotations, hmms_path)


@pytest.fixture
def indexes(run_subprocess, work_path) -> List[WFIndex]:
    index_path = work_path / "references"
    index_path.mkdir(parents=True)

    index_path = index_path / "foo"

    shutil.copytree(INDEX_PATH, index_path)

    class IndexFactory(ModelFactory):
        __model__ = Index

    return [
        WFIndex(
            IndexFactory.build(),
            index_path,
            make_mocked_coro(),
            make_mocked_coro(),
            run_subprocess,
        )
    ]


@pytest.fixture
async def sample() -> WFSample:
    faker = Faker()

    class SampleFactory(ModelFactory):
        __faker__ = faker
        __model__ = Sample

        library_type = Use(lambda: LibraryType.normal)
        paired = Use(lambda: False)

    _sample = WFSample.parse_obj(SampleFactory.build())

    return _sample


@pytest.fixture
async def reads(sample: WFSample, work_path: Path):
    reads_path = work_path / "reads"
    reads_path.mkdir(parents=True)

    shutil.copy(FASTQ_PATH, reads_path / "reads_1.fq.gz")

    return Reads(sample, {}, reads_path)


@pytest.fixture
async def analysis(
    indexes: List[Index], sample: WFSample, subtractions: List[WFSubtraction]
):
    return WFAnalysis(
        make_mocked_coro(),
        id="foo",
        created_at=arrow.utcnow().naive,
        files=[],
        index=indexes[0].index,
        job=JobNested(id="bar"),
        ready=True,
        reference=ReferenceNested(
            id="ref", data_type=ReferenceDataType.genome, name="Reference 1"
        ),
        sample=sample,
        subtractions=subtractions,
        user=UserNested(id="abc12345", handle="bob", administrator=False),
        workflow="nuvs",
    )


@pytest.fixture
async def subtractions(work_path):
    subtractions_path = work_path / "subtractions"
    subtractions_path.mkdir(parents=True)

    path_1 = work_path / "subtractions" / "arabidopsis_thaliana_1"
    path_2 = work_path / "subtractions" / "arabidopsis_thaliana_2"

    copytree(SUBTRACTION_PATH, path_1)
    copytree(SUBTRACTION_PATH, path_2)

    return [
        WFSubtraction(
            id="arabidopsis_thaliana_1",
            count=12,
            created_at=arrow.utcnow().naive,
            file=SubtractionUpload(id=12, name="arabidopsis.fa.gz"),
            files=[],
            gc=NucleotideComposition(a=0.1, t=0.2, g=0.3, c=0.4, n=0.0),
            linked_samples=[],
            path=path_1,
            name="Arabidopsis thaliana 1",
            nickname="Thalecress",
            ready=True,
            user=UserNested(administrator=False, id="bob", handle="Bob"),
        ),
        WFSubtraction(
            id="arabidopsis_thaliana_2",
            count=12,
            created_at=arrow.utcnow().naive,
            file=SubtractionUpload(id=12, name="arabidopsis.fa.gz"),
            files=[],
            gc=NucleotideComposition(a=0.1, t=0.2, g=0.3, c=0.4, n=0.0),
            linked_samples=[],
            path=path_2,
            name="Arabidopsis thaliana 2",
            nickname="Thalecress",
            ready=True,
            user=UserNested(administrator=False, id="bob", handle="Bob"),
        ),
    ]


async def test_eliminate_otus(
    caplog, indexes, reads: Reads, run_subprocess, work_path: Path
):
    caplog.set_level(logging.INFO)

    await eliminate_otus(indexes, 1, reads, run_subprocess, work_path)

    actual_path = work_path / "unmapped_otus.fq"

    with open(actual_path, "r") as f:
        actual = [line.rstrip() for line in f]
        actual = {tuple(actual[i : i + 4]) for i in range(0, len(actual), 4)}

    assert len(actual) > 0

    with open(TEST_DATA_PATH / "unmapped_otus.fq", "r") as f:
        expected = [line.rstrip() for line in f]
        expected = {tuple(expected[i : i + 4]) for i in range(0, len(expected), 4)}

    assert actual == expected


@pytest.mark.parametrize("no_subtractions", [True, False])
async def test_eliminate_subtraction(
    run_subprocess, subtractions: List[WFSubtraction], work_path, no_subtractions
):
    if no_subtractions:
        subtractions = []

    shutil.copy(TEST_DATA_PATH / "unmapped_otus.fq", work_path / "unmapped_otus.fq")
    await eliminate_subtraction(2, run_subprocess, subtractions, work_path)

    assert Path(work_path / "unmapped_subtraction.fq").is_file()

    if no_subtractions:
        with open(work_path / "unmapped_subtraction.fq") as subtracted_file:
            with open(work_path / "unmapped_otus.fq") as otu_file:
                subtracted_lines = [line.strip() for line in subtracted_file]
                otu_lines = [line.strip() for line in otu_file]
                for subtracted, otu in zip(subtracted_lines, otu_lines):
                    assert subtracted == otu


@pytest.mark.flaky(reruns=3)
@pytest.mark.parametrize("paired", [False, True], ids=["unpaired", "paired"])
async def test_reunite_pairs(paired, reads: Reads, sample: WFSample, work_path):
    if paired:
        sample.paired = paired

        unite_path = TEST_DATA_PATH / "unite.json"

        with open(unite_path, "r") as f:
            unite = json.load(f)

        for path, key in [(reads.left, "left"), (reads.right, "right")]:
            with gzip.open(path, "wt") as f:
                for line in unite[key]:
                    f.write(line + "\n")

        separate_path = work_path / "unmapped_subtractions.fq"

        with open(separate_path, "w") as f:
            for line in unite["separate"]:
                f.write(line + "\n")

        assert reads.right.exists()

    assert reads.left.exists()

    if paired:
        assert reads.right.exists()

    await reunite_pairs(reads, work_path)

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
    analysis: WFAnalysis,
    sample: Sample,
    run_subprocess,
    work_path: Path,
):
    sample.paired = paired

    proc = 1
    mem = 5

    if paired:
        for suffix in (1, 2):
            filename = f"unmapped_{suffix}.fq"

            shutil.copy(
                TEST_DATA_PATH / filename,
                work_path / filename,
            )
    else:
        shutil.copy(
            TEST_DATA_PATH / "unmapped_1.fq", work_path / "unmapped_subtractions.fq"
        )

    await assemble(analysis, mem, proc, run_subprocess, sample, work_path)

    expected_path = TEST_DATA_PATH / f"scaffolds_{'p' if paired else 'u'}.fa"
    scaffolds_path = work_path / "spades/scaffolds.fasta"
    compressed_path = work_path / "assembly.fa.gz"

    expected = {
        (record.id, record.seq) for record in SeqIO.parse(expected_path, "fasta")
    }

    # Check that scaffolds.fasta from SPAdes matches expected data.
    assert {
        (record.id, record.seq) for record in SeqIO.parse(scaffolds_path, "fasta")
    } == expected

    # Check that compressed assembly matches expected data.
    with gzip.open(compressed_path, "rt") as f:
        assert {
            (record.id, record.seq) for record in SeqIO.parse(f, "fasta")
        } == expected

    assert analysis._to_upload == [(compressed_path, "fasta")]


async def test_process_fasta(
    data_regression,
    file_regression,
    analysis: WFAnalysis,
    work_path: Path,
):
    spades_path = work_path / "spades"
    spades_path.mkdir(parents=True)

    shutil.copy(TEST_DATA_PATH / "scaffolds_u.fa", spades_path / "scaffolds.fasta")

    results = dict()

    await process_fasta(analysis, 2, results, work_path)

    data_regression.check(results)

    with open(work_path / "orfs.fa") as f:
        file_regression.check(f.read())


async def test_vfam(
    data_regression, analysis: WFAnalysis, hmms: HMMs, run_subprocess, work_path: Path
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
