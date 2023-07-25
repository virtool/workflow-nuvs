import asyncio
import collections
import shlex
import shutil
from pathlib import Path
from typing import List

import aiofiles
import rust
from virtool_core.bio import (
    find_orfs,
    read_fasta,
)
from virtool_core.utils import compress_file
from virtool_workflow import hooks
from virtool_workflow import step
from virtool_workflow.analysis.hmms import HMMs
from virtool_workflow.analysis.reads import Reads
from virtool_workflow.data_model.analysis import WFAnalysis
from virtool_workflow.data_model.indexes import WFIndex
from virtool_workflow.data_model.subtractions import WFSubtraction


@hooks.on_failure
async def delete_analysis_document(analysis_provider):
    await analysis_provider.delete()


@hooks.on_success
async def upload_result(analysis_provider, results):
    await analysis_provider.upload_result(results)


@step(name="Eliminate OTUs")
async def eliminate_otus(
    indexes: List[WFIndex], proc: int, reads: Reads, run_subprocess, work_path: Path
):
    """
    Map sample reads to reference OTUs and discard.

    Bowtie2 is set to use the search parameter ``--very-fast-local`` and retain
    unaligned reads to the FASTQ file ``unmapped_subtraction.fq``.

    """
    index_path = str(indexes[0].bowtie_path)

    paths = [str(reads.left)]

    if reads.sample.paired:
        paths.append(str(reads.right))

    command = [
        "bowtie2",
        "-p",
        str(proc),
        "-k",
        str(1),
        "--very-fast-local",
        "-x",
        index_path,
        "--un",
        str(work_path / "unmapped_otus.fq"),
        "-U",
        ",".join(paths),
    ]

    await run_subprocess(command)


@step
async def eliminate_subtraction(
    proc: int, run_subprocess, subtractions: List[WFSubtraction], work_path: Path
):
    """
    Map remaining reads to the subtraction and discard.

    Reads that were not mapped to the reference OTUs in the previous step
    (`unmapped_otus.fq`) are mapped against the subtraction. Reads with no
    alignment against the subtraction (`unmapped_subtractions.fq`) are carried
    forward into the next step.

    Bowtie2 is set to use the search parameter ``--very-fast-local`` and retain
    unaligned reads to the FASTQ file ``unmapped_subtraction.fq``. Providing the `--un`
    option to Bowtie2 writes any unmapped reads to the path provided with the
    option.

    """

    if len(subtractions) == 0:
        await asyncio.to_thread(
            shutil.copyfile,
            work_path / "unmapped_otus.fq",
            work_path / "unmapped_subtraction.fq",
        )

    await asyncio.to_thread(
        shutil.copyfile, work_path / "unmapped_otus.fq", work_path / "working_otus.fq"
    )

    for subtraction in subtractions:
        command = [
            "bowtie2",
            "--very-fast-local",
            "-k",
            str(1),
            "-p",
            str(proc),
            "-x",
            shlex.quote(str(subtraction.bowtie2_index_path)),
            "--un",
            str(work_path / "unmapped_subtraction.fq"),
            "-U",
            str(work_path / "working_otus.fq"),
        ]

        await run_subprocess(command)

        await asyncio.to_thread(
            shutil.copyfile,
            work_path / "unmapped_subtraction.fq",
            work_path / "working_otus.fq",
        )


@step
async def reunite_pairs(reads: Reads, work_path: Path):
    """
    Reunite paired reads after elimination.
    """
    if reads.sample.paired:
        rs_reads = rust.Reads(
            reads.sample.paired,
            str(reads.left.absolute()),
            str(reads.right.absolute()),
            str(work_path / "unmapped_subtractions.fq"),
        )

        rust.reunite_pairs(rs_reads, str(work_path))


@step
async def assemble(
    analysis: WFAnalysis,
    mem: int,
    proc: int,
    run_subprocess,
    sample,
    work_path: Path,
):
    """Assemble reads using SPAdes."""
    command = ["spades.py", "-t", str(proc), "-m", str(mem)]

    if sample.paired:
        command += [
            "-1",
            str(work_path / "unmapped_1.fq"),
            "-2",
            str(work_path / "unmapped_2.fq"),
        ]
    else:
        command += [
            "-s",
            str(work_path / "unmapped_subtractions.fq"),
        ]

    k = "21,33,55,75"

    if sample.library_type == "srna":
        k = "17,21,23"

    spades_path = work_path / "spades"

    command += ["-o", str(spades_path), "-k", k]

    await run_subprocess(command)

    compressed_assembly_path = work_path / "assembly.fa.gz"

    await asyncio.to_thread(
        compress_file,
        spades_path / "scaffolds.fasta",
        compressed_assembly_path,
        processes=proc,
    )

    analysis.upload(compressed_assembly_path, "fasta")


@step
async def process_fasta(
    analysis: WFAnalysis,
    proc: int,
    results: dict,
    work_path: Path,
):
    """
    Find ORFs in the assembled contigs.

    Only ORFs that are 100+ amino acids long are recorded. Contigs with no acceptable
    ORFs are discarded.

    """
    assembly_path = work_path / "spades/scaffolds.fasta"

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
            orf["hits"] = list()

        # Make an entry for the nucleotide sequence containing a unique integer index,
        # the sequence itself, and all ORFs in the sequence.
        sequences.append({"index": len(sequences), "sequence": sequence, "orfs": orfs})

    # Write the ORFs to a FASTA file so that they can be analyzed using HMMER and vFAM.
    orfs_path = work_path / "orfs.fa"

    async with aiofiles.open(orfs_path, "w") as f:
        for entry in sequences:
            for orf in entry["orfs"]:
                await f.write(
                    f">sequence_{entry['index']}.{orf['index']}\n{orf['pro']}\n"
                )

    compressed_orfs_path = Path(f"{orfs_path}.gz")

    await asyncio.to_thread(
        compress_file, orfs_path, compressed_orfs_path, processes=proc
    )

    analysis.upload(compressed_orfs_path, "fasta")

    results["hits"] = sequences


@step(name="VFam")
async def vfam(
    analysis: WFAnalysis,
    hmms: HMMs,
    proc: int,
    results: dict,
    run_subprocess,
    work_path: Path,
):
    """
    Search for viral motifs in ORF translations.

    ORF translations are generated by :meth:`.process_fasta`. Viral motifs are found
    using ``hmmscan`` to search through ``candidates.fa`` using the profile HMMs in
    ``data_path/hmm/vFam.hmm``.

    Saves two files:

    - ``hmm.tsv`` contains the raw output of `hmmer`
    - ``hits.tsv`` contains the `hmmer` results formatted and annotated with the
      annotations from the Virtool HMM database collection

    """
    tsv_path = work_path / "hmm.tsv"

    command = [
        "hmmscan",
        "--tblout",
        str(tsv_path),
        "--noali",
        "--cpu",
        str(proc - 1),
        str(hmms.path / "profiles.hmm"),
        str(work_path / "orfs.fa"),
    ]

    await run_subprocess(command)

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
                    }
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

    analysis.upload(tsv_path, "tsv")
