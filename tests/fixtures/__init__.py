from pathlib import Path

import pytest

from tests.fixtures.unite import unite

__all__ = ["unite"]


@pytest.fixture()
def example_path():
    return Path(__file__).parent.parent.parent / "example"


@pytest.fixture()
def virtool_workflow_example_path(example_path: Path):
    return example_path


@pytest.fixture()
def work_path(tmpdir):
    path = Path(tmpdir) / "work"
    path.mkdir(parents=True)

    return path
