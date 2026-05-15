"""
Pipeline-side helpers for virusHanter2.

These are the small utilities that are called directly from the Snakemake
rule files: sample discovery, FASTQ/FASTA wrangling, Kraken2 post-processing,
and the BLASTN driver. All HTML-report and plotting helpers that used to live
here have moved into the reportHanter package (`reporthanter` on PyPI), which
the report rules invoke either by CLI or by library import.
"""

from pathlib import Path
import os
import re
import subprocess

import numpy as np
import pandas as pd
import pyfastx


def read_file_as_blob(file_path: str) -> str:
    """Return the file content as a lowercase hex string.

    Used by `aggregate_run_information` to embed per-sample HTML reports as
    a column in the run summary CSV (matches the original virusHanter
    behavior).
    """
    with open(file_path, "rb") as f:
        return f.read().hex()


def common_suffix(folder: str) -> str:
    """Find the longest filename suffix shared by every sequencing file in
    `folder`. Used to derive the read1/read2 file extension at workflow
    setup time.
    """
    samples = sorted(
        file.name
        for file in Path(folder).iterdir()
        if re.search(r"\.(fq|fastq|fa|fasta|fna)$", file.name)
    )
    if not samples:
        return ""

    test_sample = samples[0]
    suffix = ""
    for i in range(1, len(test_sample) + 1):
        index = -i
        if any(sample[index] != test_sample[index] for sample in samples):
            break
        suffix += test_sample[index]
    return suffix[::-1]


def paired_reads(folder: str) -> list:
    """Return the common-prefix sample names for paired-end FASTQ files in
    `folder`. Assumes files are sorted by name and pair up two-by-two.
    """
    def common_prefix(a: str, b: str) -> str:
        out = ""
        for ca, cb in zip(a, b):
            if ca != cb:
                break
            out += ca
        return out

    samples = sorted(
        x.stem
        for x in Path(folder).iterdir()
        if re.search(r"\.(fq|fastq|fa|fasta|fna)$", x.name)
    )

    prefixes = []
    for i in range(0, len(samples), 2):
        prefixes.append(common_prefix(samples[i], samples[i + 1]))
    return prefixes


def kaiju_db_files(kaiju_db: str) -> tuple:
    """Locate the `.fmi`, `names.dmp`, and `nodes.dmp` files inside a Kaiju
    database directory.

    Returns placeholder paths under the given directory when the directory
    itself does not yet exist. This lets workflow construction (including
    `snakemake --lint` and `snakemake -n`) proceed even before the database
    has been materialized; the rule that consumes these paths will still
    fail loudly at run time if the files are not present.
    """
    db_path = Path(kaiju_db)
    if not db_path.is_dir():
        return (
            db_path / "kaiju_db.fmi",
            db_path / "names.dmp",
            db_path / "nodes.dmp",
        )

    files = [x for x in db_path.iterdir() if x.is_file()]
    fmi = next((x for x in files if x.suffix == ".fmi"), db_path / "kaiju_db.fmi")
    names = next((x for x in files if x.name == "names.dmp"), db_path / "names.dmp")
    nodes = next((x for x in files if x.name == "nodes.dmp"), db_path / "nodes.dmp")
    return fmi, names, nodes


def fastx_file_to_df(fastx_file: str) -> pd.DataFrame:
    """Read a FASTA/FASTQ file into a DataFrame sorted by sequence length.

    pyfastx 2.x yields tuples from `Fastx`; the first two fields are
    (name, sequence). Earlier versions returned attribute-bearing objects;
    callers in the original virusHanter targeted that older API.
    """
    fastx = pyfastx.Fastx(fastx_file)
    rows = [(record[0], record[1]) for record in fastx]
    if not rows:
        return pd.DataFrame(columns=["name", "sequence", "read_len"])

    names, seqs = zip(*rows)
    return (
        pd.DataFrame({"name": list(names), "sequence": list(seqs)})
        .assign(read_len=lambda x: x.sequence.str.len())
        .sort_values("read_len", ascending=False)
    )


def wrangle_kraken(kraken_file: str) -> pd.DataFrame:
    """Parse a Kraken2 report TSV into a DataFrame with an explicit `domain`
    column carried down from the nearest D/U/R parent row.
    """
    kraken = (
        pd.read_csv(
            kraken_file,
            sep="\t",
            header=None,
            names=["percent", "count_clades", "count", "tax_lvl", "taxonomy_id", "name"],
        )
        .assign(name=lambda x: x.name.str.strip())
        .assign(
            domain=lambda x: np.select(
                [x.tax_lvl.isin(["D", "U", "R"])],
                [x.name],
                default=pd.NA,
            )
        )
    )
    kraken["domain"] = kraken["domain"].ffill()
    return kraken


def run_blastn(contigs_csv: str, db: str, temp_file: str, threads: int) -> pd.DataFrame:
    """Run blastn (megablast) on each row of a contigs CSV one at a time and
    return the input table joined with the best hit per contig.
    """
    os.environ["BLASTDB"] = db
    df = pd.read_csv(contigs_csv)
    if df.empty:
        return df

    matches = []
    for contig in df.itertuples():
        with open(temp_file, "w") as f:
            f.write(f">{contig.name}\n{contig.sequence}\n")
        command = [
            "blastn", "-num_threads", str(threads), "-task", "megablast",
            "-query", temp_file, "-db", db, "-max_target_seqs", "1",
            "-outfmt", "6 stitle sacc pident slen",
        ]
        # capture_output so a non-zero blastn includes stderr in the raised
        # CalledProcessError rather than printing to the main process stdout.
        result = subprocess.run(
            command,
            capture_output=True,
            text=True,
            check=True,
        )
        matches.append(result.stdout.strip())

    df = df.assign(matches=matches).loc[lambda x: x.matches != ""]
    if df.empty:
        return df

    df[["match_name", "accession", "percent_identity", "sequence_len"]] = (
        df.matches.str.split("\t", expand=True).iloc[:, :4]
    )
    df = df.assign(sequence_len=lambda x: x.sequence_len.str.split("\n").str[0])
    return df
