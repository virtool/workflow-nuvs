import collections
import shlex
from pathlib import Path
from typing import List

import aiofiles
from virtool_core.bio import read_fastq_headers, read_fastq_from_path, find_orfs, read_fasta
from virtool_workflow import step
from virtool_workflow.analysis.analysis import Analysis
from virtool_workflow.analysis.hmms import HMMs
from virtool_workflow.analysis.indexes import Index
from virtool_workflow.analysis.reads import Reads
from virtool_workflow.execution.run_in_executor import FunctionExecutor


@step
async def eliminate_otus(indexes: List[Index], proc: int, reads: Reads, run_subprocess, work_path: Path):
    index_path = str(indexes[0].bowtie_path)

    command = [
        "bowtie2",
        "-p", str(proc),
        "-k", str(1),
        "--very-fast-local",
        "-x", index_path,
        "--un", str(work_path / "unmapped_otus.fq"),
        "-U", ",".join([str(path) for path in reads.paths])
    ]

    await run_subprocess(command)


@step
async def eliminate_subtraction(proc: int, run_subprocess, subtractions, work_path: Path):
    """
    Map unaligned reads the `eliminate_otus` step sample's subtraction host using ``bowtie2``.

    Bowtie2 is set to use the search parameter ``--very-fast-local`` and retain unaligned reads to the FASTQ file
    ``unmapped_host.fq``.

    """
    command = [
        "bowtie2",
        "--very-fast-local",
        "-k", str(1),
        "-p", str(proc),
        "-x", shlex.quote(subtractions[0].path),
        "--un", str(work_path / "unmapped_hosts.fq"),
        "-U", str(work_path / "unmapped_otus.fq")
    ]

    await run_subprocess(command)


@step
async def reunite_pairs(reads, sample, work_path):
    if sample.paired:
        unmapped_path = work_path / "unmapped_hosts.fq"
        headers = await read_fastq_headers(unmapped_path)

        unmapped_roots = {h.split(" ")[0] for h in headers}

        async with aiofiles.open(work_path / "unmapped_1.fq", "w") as f:
            async for header, seq, quality in read_fastq_from_path(reads.paths[0]):
                if header.split(" ")[0] in unmapped_roots:
                    await f.write("\n".join([header, seq, "+", quality]) + "\n")

        async with aiofiles.open(work_path / "unmapped_2.fq", "w") as f:
            async for header, seq, quality in read_fastq_from_path(reads.paths[1]):
                if header.split(" ")[0] in unmapped_roots:
                    await f.write("\n".join([header, seq, "+", quality]) + "\n")


@step
async def assemble(analysis, mem, proc, run_subprocess, run_in_executor, sample, work_path):
    """
    Call ``spades.py`` to assemble contigs from ``unmapped_hosts.fq``. Passes ``21,33,55,75`` for the ``-k``
    argument.

    """
    command = [
        "spades.py",
        "-t", str(proc),
        "-m", str(mem)
    ]

    if sample.paired:
        command += [
            "-1", str(work_path / "unmapped_1.fq"),
            "-2", str(work_path / "unmapped_2.fq"),
        ]
    else:
        command += [
            "-s", str(work_path / "unmapped_hosts.fq"),
        ]

    k = "21,33,55,75"

    if sample.library_type == "srna":
        k = "17,21,23"

    spades_path = work_path / "spades"

    command += [
        "-o", str(spades_path),
        "-k", k
    ]

    await run_subprocess(command)

    await analysis.upload_file(
        spades_path / "scaffolds.fasta",
    )


@step
async def process_fasta(analysis: Analysis, results: dict, run_in_executor: FunctionExecutor, work_path: Path):
    """
    Finds ORFs in the contigs assembled by :meth:`.assemble`. Only ORFs that are 100+ amino acids long are recorded.
    Contigs with no acceptable ORFs are discarded.

    """
    assembly_path = work_path / "assembly.fa"

    assembly = await run_in_executor(
        read_fasta,
        assembly_path
    )

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

        # Make an entry for the nucleotide sequence containing a unique integer index, the sequence itself, and
        # all ORFs in the sequence.
        sequences.append({
            "index": len(sequences),
            "sequence": sequence,
            "orfs": orfs
        })

    # Write the ORFs to a FASTA file so that they can be analyzed using HMMER and vFAM.
    orfs_path = work_path / "orfs.fa"

    async with aiofiles.open(orfs_path, "w") as f:
        for entry in sequences:
            for orf in entry["orfs"]:
                await f.write(f">sequence_{entry['index']}.{orf['index']}\n{orf['pro']}\n")

    await analysis.upload_file(
        orfs_path
    )

    results["sequences"] = sequences


@step
async def vfam(hmms: HMMs, proc: int, results: dict, run_subprocess, work_path: Path):
    """
    Searches for viral motifs in ORF translations generated by :meth:`.process_fasta`. Calls ``hmmscan`` and
    searches against ``candidates.fa`` using the profile HMMs in ``data_path/hmm/vFam.hmm``.

    Saves two files:

    - ``hmm.tsv`` contains the raw output of `hmmer`
    - ``hits.tsv`` contains the `hmmer` results formatted and annotated with the annotations from the Virtool HMM
      database collection

    """
    # The path to output the hmmer results to.
    tsv_path = work_path / "hmm.tsv"

    command = [
        "hmmscan",
        "--tblout", tsv_path,
        "--noali",
        "--cpu", str(proc - 1),
        hmms.path,
        work_path / "orfs.fa"
    ]

    await run_subprocess(command)

    hits = collections.defaultdict(lambda: collections.defaultdict(list))

    # Go through the raw HMMER results and annotate the HMM hits with data from the database.
    async with aiofiles.open(tsv_path, "r") as f:
        async for line in f:
            if line.startswith("vFam"):
                line = line.split()

                cluster_id = int(line[0].split("_")[1])

                annotation_id = hmms.cluster_annotation_map[cluster_id]

                # Expecting sequence_0.0
                sequence_index, orf_index = (int(x) for x in line[2].split("_")[1].split("."))

                hits[sequence_index][orf_index].append({
                    "hit": annotation_id,
                    "full_e": float(line[4]),
                    "full_score": float(line[5]),
                    "full_bias": float(line[6]),
                    "best_e": float(line[7]),
                    "best_bias": float(line[8]),
                    "best_score": float(line[9])
                })

    sequences = results["sequences"]

    for sequence_index in hits:
        for orf_index in hits[sequence_index]:
            sequences[sequence_index]["orfs"][orf_index]["hits"] = hits[sequence_index][orf_index]

        sequence = sequences[sequence_index]

        if all(len(o["hits"]) == 0 for o in sequence["orfs"]):
            sequences.remove(sequence)
