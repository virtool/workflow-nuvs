import asyncio
import collections
import os
import shlex
import shutil
from pathlib import Path

import aiofiles
from pyfixtures import fixture
from structlog import get_logger
from virtool.bio import (
    find_orfs,
    read_fasta,
)
from virtool.models.enums import LibraryType
from virtool.utils import compress_file, decompress_file
from virtool_workflow import RunSubprocess, hooks, step
from virtool_workflow.analysis.skewer import (
    SkewerConfiguration,
    SkewerMode,
    SkewerRunner,
)
from virtool_workflow.analysis.trimming import calculate_trimming_min_length
from virtool_workflow.analysis.utils import ReadPaths
from virtool_workflow.data.analyses import WFAnalysis
from virtool_workflow.data.hmms import WFHMMs
from virtool_workflow.data.indexes import WFIndex
from virtool_workflow.data.samples import WFSample
from virtool_workflow.data.subtractions import WFSubtraction

from utils import filter_reads_by_headers, read_fastq_headers


@hooks.on_failure
async def delete_analysis(analysis: WFAnalysis):
    await analysis.delete()


@fixture
async def trimmed_path(work_path: Path) -> Path:
    trimmed_path = work_path / "trimmed"
    trimmed_path.mkdir(exist_ok=True)

    return trimmed_path


@fixture
async def trimmed_read_paths(sample: WFSample, trimmed_path: Path) -> ReadPaths:
    if sample.paired:
        return (
            trimmed_path / "reads_1.fq.gz",
            trimmed_path / "reads_2.fq.gz",
        )

    return (trimmed_path / "reads_1.fq.gz",)


@step()
async def trim_reads(
    proc: int,
    sample: WFSample,
    skewer: SkewerRunner,
    work_path: Path,
):
    """Trim reads using Skewer"""
    trimmed_path = work_path / "trimmed"
    await asyncio.to_thread(trimmed_path.mkdir, parents=True)

    await skewer(
        SkewerConfiguration(
            min_length=calculate_trimming_min_length(sample),
            mode=SkewerMode.PAIRED_END if sample.paired else SkewerMode.SINGLE_END,
            number_of_processes=proc,
        ),
        sample.read_paths,
        output_path=trimmed_path,
    )


@step(name="Eliminate OTUs")
async def eliminate_otus(
    index: WFIndex,
    proc: int,
    run_subprocess: RunSubprocess,
    trimmed_read_paths: ReadPaths,
    work_path: Path,
):
    """Map sample reads to reference OTUs and discard.

    Bowtie2 is set to use the search parameter ``--very-fast-local`` and retain
    unaligned reads to the FASTQ file ``unmapped_subtraction.fq``.

    """
    command = [
        "bowtie2",
        "-p",
        proc,
        "-k",
        1,
        "--very-fast-local",
        "-x",
        index.bowtie_path,
        "--un",
        work_path / "unmapped_otus.fq",
        "-U",
        *trimmed_read_paths,
    ]

    await run_subprocess(command)


@step
async def eliminate_subtraction(
    proc: int,
    run_subprocess: RunSubprocess,
    subtractions: list[WFSubtraction],
    work_path: Path,
):
    """Map remaining reads to the subtraction and discard.

    Reads that were not mapped to the reference OTUs in the previous step
    (`unmapped_otus.fq`) are mapped against the subtraction. Reads with no
    alignment against the subtraction (`unmapped_subtractions.fq`) are carried
    forward into the next step.

    Bowtie2 is set to use the search parameter ``--very-fast-local`` and retain
    unaligned reads to the FASTQ file ``unmapped_subtraction.fq``. Providing the `--un`
    option to Bowtie2 writes any unmapped reads to the path provided with the
    option.

    """
    if subtractions:
        await asyncio.to_thread(
            shutil.copyfile,
            work_path / "unmapped_otus.fq",
            work_path / "working_otus.fq",
        )

        for subtraction in subtractions:
            await run_subprocess(
                [
                    "bowtie2",
                    "--very-fast-local",
                    "-k",
                    1,
                    "-p",
                    proc,
                    "-x",
                    shlex.quote(str(subtraction.bowtie2_index_path)),
                    "--un",
                    work_path / "unmapped_subtractions.fq",
                    "-U",
                    work_path / "working_otus.fq",
                ],
            )

            await asyncio.to_thread(
                shutil.copyfile,
                work_path / "unmapped_subtractions.fq",
                work_path / "working_otus.fq",
            )

        await asyncio.to_thread(
            os.rename,
            work_path / "working_otus.fq",
            work_path / "unmapped_subtractions.fq",
        )

    else:
        await asyncio.to_thread(
            shutil.copyfile,
            work_path / "unmapped_otus.fq",
            work_path / "unmapped_subtractions.fq",
        )


@step
async def reunite_pairs(
    proc: int,
    sample: WFSample,
    trimmed_read_paths: ReadPaths,
    work_path: Path,
):
    """Reunite paired reads after elimination."""
    if sample.paired:
        headers = await asyncio.to_thread(
            read_fastq_headers,
            work_path / "unmapped_subtractions.fq",
        )

        for path in trimmed_read_paths:
            await asyncio.to_thread(
                decompress_file,
                path,
                path.with_suffix(".fq"),
                proc,
            )

        path_1, path_2 = trimmed_read_paths

        await asyncio.to_thread(
            filter_reads_by_headers,
            headers,
            (
                work_path / "unmapped_1.fq",
                work_path / "unmapped_2.fq",
            ),
            (path_1.with_suffix(".fq"), path_2.with_suffix(".fq")),
        )


@step
async def assemble(
    analysis: WFAnalysis,
    mem: int,
    proc: int,
    run_subprocess: RunSubprocess,
    sample: WFSample,
    work_path: Path,
):
    """Assemble reads using SPAdes."""
    spades_path = work_path / "spades"

    k = "21,33,55,75"

    if sample.library_type == LibraryType.srna:
        k = "17,21,23"

    command = [
        "spades.py",
        "-t",
        proc,
        "-m",
        mem,
        "-k",
        k,
        "-o",
        spades_path,
    ]

    logger = get_logger("spades")

    if sample.paired:
        command += [
            "-1",
            work_path / "unmapped_1.fq",
            "-2",
            work_path / "unmapped_2.fq",
        ]
    else:
        command += [
            "-s",
            work_path / "unmapped_subtractions.fq",
        ]

    async def handler(line):
        logger.info("stdout", line=line.decode().strip())

    await run_subprocess(command, stdout_handler=handler)

    compressed_assembly_path = work_path / "assembly.fa.gz"

    await asyncio.to_thread(
        compress_file,
        spades_path / "scaffolds.fasta",
        compressed_assembly_path,
        processes=proc,
    )

    await analysis.upload_file(compressed_assembly_path, "fasta")


@step
async def process_assembly(
    analysis: WFAnalysis,
    proc: int,
    results: dict,
    work_path: Path,
):
    """Find ORFs in the assembled contigs.

    Only ORFs that are 100+ amino acids long are recorded. Contigs with no acceptable
    ORFs are discarded.

    """
    assembly_path = work_path / "spades/scaffolds.fa"

    await asyncio.to_thread(
        os.rename,
        work_path / "spades/scaffolds.fasta",
        assembly_path,
    )

    assembly = await asyncio.to_thread(read_fasta, assembly_path)

    sequences = list()

    for _, sequence in assembly:
        sequence_length = len(sequence)

        # Don't consider the sequence if it is shorter than 300 bp.
        if sequence_length < 300:
            continue

        orfs = find_orfs(sequence)

        # Don't consider the sequence if it has no ORFs.
        if len(orfs) == 0:
            continue

        # Add an index field to each orf dict.
        orfs = [dict(o, index=i) for i, o in enumerate(orfs)]

        for orf in orfs:
            orf.pop("nuc")
            orf["hits"] = []

        # Make an entry for the nucleotide sequence containing a unique integer index,
        # the sequence itself, and all ORFs in the sequence.
        sequences.append({"index": len(sequences), "sequence": sequence, "orfs": orfs})

    # Write the ORFs to a FASTA file so that they can be analyzed using HMMER and vFAM.
    orfs_path = work_path / "orfs.fa"

    async with aiofiles.open(orfs_path, "w") as f:
        for entry in sequences:
            for orf in entry["orfs"]:
                await f.write(
                    f">sequence_{entry['index']}.{orf['index']}\n{orf['pro']}\n",
                )

    compressed_orfs_path = Path(f"{orfs_path}.gz")

    await asyncio.to_thread(
        compress_file,
        orfs_path,
        compressed_orfs_path,
        processes=proc,
    )

    await analysis.upload_file(compressed_orfs_path, "fasta")

    results["hits"] = sequences


@step(name="VFam")
async def vfam(
    analysis: WFAnalysis,
    hmms: WFHMMs,
    proc: int,
    results: dict,
    run_subprocess: RunSubprocess,
    work_path: Path,
):
    """Search for viral motifs in ORF translations.

    ORF translations are generated by :meth:`.process_fasta`. Viral motifs are found
    using ``hmmscan`` to search through ``candidates.fa`` using the profile HMMs in
    ``data_path/hmm/vFam.hmm``.

    Saves two files:

    - ``hmm.tsv`` contains the raw output of `hmmer`
    - ``hits.tsv`` contains the `hmmer` results formatted and annotated with the
      annotations from the Virtool HMM database collection

    """
    tsv_path = work_path / "hmm.tsv"

    await run_subprocess(
        [
            "hmmscan",
            "--tblout",
            tsv_path,
            "--noali",
            "--cpu",
            proc - 1,
            hmms.path / "profiles.hmm",
            work_path / "orfs.fa",
        ],
    )

    hmmer_hits = collections.defaultdict(lambda: collections.defaultdict(list))

    # Go through the raw HMMER results and annotate the HMM hits with data from the
    # database.
    async with aiofiles.open(tsv_path, "r") as f:
        async for line in f:
            if line.startswith("vFam"):
                line = line.split()

                cluster_id = int(line[0].split("_")[1])

                annotation_id = hmms.cluster_annotation_map[cluster_id]

                # Expecting sequence_0.0
                sequence_index, orf_index = (
                    int(x) for x in line[2].split("_")[1].split(".")
                )

                hmmer_hits[sequence_index][orf_index].append(
                    {
                        "hit": annotation_id,
                        "full_e": float(line[4]),
                        "full_score": float(line[5]),
                        "full_bias": float(line[6]),
                        "best_e": float(line[7]),
                        "best_bias": float(line[8]),
                        "best_score": float(line[9]),
                    },
                )

    hits = results["hits"]

    for sequence_index in hmmer_hits:
        for orf_index in hmmer_hits[sequence_index]:
            hits[sequence_index]["orfs"][orf_index]["hits"] = hmmer_hits[
                sequence_index
            ][orf_index]

        sequence = results["hits"][sequence_index]

        if all(len(orf["hits"]) == 0 for orf in sequence["orfs"]):
            hits.remove(sequence)

    await analysis.upload_file(tsv_path, "tsv")
    await analysis.upload_result(results)
