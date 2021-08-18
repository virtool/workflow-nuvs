from pathlib import Path

from setuptools import setup, find_packages

AUTHORS = ["Ian Boyes"]

CLASSIFIERS = [
    "Programming Language :: Python :: 3.9",
]

PACKAGES = find_packages(exclude="tests")

INSTALL_REQUIRES = [
    "aiofiles>=0.6.0",
    "virtool-core>=0.3.0",
    "virtool-workflow==0.6.0"
]

setup(
    name="vt-workflow-nuvs",
    version="0.1.0",
    description="A workflow for detecting viruses in high-throughput sequencing (HTS) libraries.",
    long_description=Path("README.md").read_text(),
    long_description_content_type="text/markdown",
    url="https://github.com/virtool/workflow-pathoscope",
    author=", ".join(AUTHORS),
    license="MIT",
    platforms="linux",
    packages=PACKAGES,
    install_requires=INSTALL_REQUIRES,
    python_requires=">=3.9",
)