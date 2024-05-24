import itertools
from pathlib import Path

from Bio import SeqIO
from virtool_workflow.analysis.utils import ReadPaths


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

    for path, handle in zip(read_paths, out_handles):
        records = (
            record for record in SeqIO.parse(path, "fastq") if record.id in headers
        )

        while True:
            batch = list(itertools.islice(records, 1000))

            if not batch:
                break

            SeqIO.write(batch, handle, "fastq")
