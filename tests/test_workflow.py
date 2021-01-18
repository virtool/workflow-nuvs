import filecmp
import json
import operator
import shutil
import sys
from shutil import copytree
from typing import List

import pytest
from pathlib import Path

from virtool_workflow.abc import AbstractFileUploader
from virtool_workflow.analysis.analysis import Analysis
from virtool_workflow.analysis.indexes import Index, Reference
from virtool_workflow.analysis.reads import Reads

import virtool_workflow.execution.run_subprocess
from virtool_workflow.analysis.subtractions.subtraction import Subtraction

from workflow import eliminate_otus, eliminate_subtraction, reunite_pairs, process_fasta

TEST_FILES_PATH = Path(sys.path[0]) / "tests"
FASTQ_PATH = TEST_FILES_PATH / "test.fq"
INDEX_PATH = TEST_FILES_PATH / "index"
NUVS_PATH = TEST_FILES_PATH / "nuvs"
SUBTRACTION_PATH = TEST_FILES_PATH / "subtraction"


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
def indexes(run_in_executor, run_subprocess, work_path) -> List[Index]:
    index_path = work_path / "references"
    index_path.mkdir(parents=True)

    index_path = index_path / "foo"

    shutil.copytree(INDEX_PATH, index_path)

    reference = Reference("bar", "genome", "Reference 1")

    index = Index("foo", index_path, reference, run_in_executor, run_subprocess)

    return [index]


@pytest.fixture
async def reads(work_path):
    reads_path = work_path / "reads"
    reads_path.mkdir(parents=True)

    shutil.copy(FASTQ_PATH, reads_path)

    return Reads(False, 75, 75, 10000, (reads_path / "test.fq",))


@pytest.fixture
async def subtractions(work_path):
    subtractions_path = work_path / "subtractions"
    subtractions_path.mkdir(parents=True)

    subtraction_path = work_path / "subtractions" / "subtraction"

    copytree(SUBTRACTION_PATH, subtraction_path)

    subtraction = Subtraction(
        "Arabidopsis thaliana",
        "Thalecress",
        subtraction_path,
        subtraction_path / "subtraction.fa",
        subtraction_path / "subtraction" / "subtraction",
        12,
        {"a": 0.1, "t": 0.2, "g": 0.3, "c": 0.4},
    )

    return [subtraction]


async def test_eliminate_otus(indexes, reads, run_subprocess, work_path):
    await eliminate_otus(indexes, 2, reads, run_subprocess, work_path)

    actual_path = work_path / "unmapped_otus.fq"

    with open(actual_path, "r") as f:
        actual = [line.rstrip() for line in f]
        actual = {tuple(actual[i : i + 4]) for i in range(0, len(actual), 4)}

    expected_path = NUVS_PATH / "unmapped_otus.fq"

    with open(expected_path, "r") as f:
        expected = [line.rstrip() for line in f]
        expected = {tuple(expected[i : i + 4]) for i in range(0, len(expected), 4)}

    assert actual == expected


async def test_eliminate_subtraction(run_subprocess, subtractions, work_path):
    shutil.copy(NUVS_PATH / "unmapped_otus.fq", work_path / "unmapped_otus.fq")

    await eliminate_subtraction(2, run_subprocess, subtractions, work_path)
