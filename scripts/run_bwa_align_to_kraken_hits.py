"""Snakemake script: build the per-sample viral reference set and
align host-removed reads to it.

Aggregates the top viral taxids from KRAKEN, KAIJU and BLAST,
walks each taxid through the parquet (with optional genus
walk-up via the NCBI taxdump), writes a multi-FASTA of the
selected references and the matching ``virus_names`` sidecar,
then runs BWA to produce the alignment BAM consumed by mosdepth
and the per-virus metrics rule.

Lives in ``scripts/`` so the ``conda:`` env declared on the
rule (``envs/bwa.yaml``) is honoured at execution time.
"""

import shlex
import shutil
import subprocess
import sys
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from scripts.functions import (  # noqa: E402
    _format_aliases,
    _load_taxdump_for_species_walkup,
    find_genus_taxid,
    find_species_taxid,
    parse_nodes_dmp,
    parquet_accession_to_taxid,
)

snakemake = snakemake  # type: ignore[name-defined]

params = snakemake.params
input_ = snakemake.input
output = snakemake.output
threads = snakemake.threads
log_path = snakemake.log[0] if snakemake.log else "/dev/null"


# Resolve bash at import time so non-standard installations are handled
# transparently. Falls back to /bin/bash if bash is not found on PATH.
_BASH = shutil.which("bash") or "/bin/bash"


def _shell(cmd: str) -> None:
    subprocess.run(cmd, shell=True, check=True, executable=_BASH)


virus_db_df = pd.read_parquet(params.virus_db)
parquet_tax_ids = set(virus_db_df["tax_id"].dropna().astype(int).tolist())
acc_to_tax = parquet_accession_to_taxid(virus_db_df)

genus_to_rep_taxid: dict[int, int] = {}
if "genus_taxid" in virus_db_df.columns:
    sub = virus_db_df.loc[virus_db_df["genus_taxid"] > 0].copy()
    if not sub.empty:
        sub = sub.assign(_seqlen=sub["sequence"].str.len())
        sub = sub.sort_values(
            ["genus_taxid", "_seqlen"],
            ascending=[True, False],
            kind="mergesort",
        )
        for r in sub.drop_duplicates(subset=["genus_taxid"], keep="first").itertuples():
            genus_to_rep_taxid[int(r.genus_taxid)] = int(r.tax_id)

rank_filter: set[str] = set(params.coverage_rank_filter or [])
taxdump_path = params.taxdump_nodes
nodes: dict[int, tuple[int, str]] = {}
sci_name_by_tid: dict[int, str] = {}
alias_by_tid: dict[int, list[str]] = {}
if taxdump_path and Path(taxdump_path).is_file():
    nodes = parse_nodes_dmp(Path(taxdump_path))
    # ``names.dmp`` is published next to ``nodes.dmp`` by the refresh
    # workflow; load the scientific name plus the alias categories
    # (acronym / common name / equivalent name / ...) so the
    # virus_names sidecar can carry both the ICTV binomial and the
    # legacy names scientists still recognise.
    names_path = Path(taxdump_path).parent / "names.dmp"
    parsed = _load_taxdump_for_species_walkup(str(taxdump_path), str(names_path))
    if parsed is not None:
        _node_info, sci_name_by_tid, alias_by_tid = parsed
elif rank_filter or params.coverage_genus_walkup:
    print(
        "[bwa_align_to_kraken_hits] TAXDUMP_NODES not set or missing; "
        "rank filter and genus walk-up are disabled."
    )
    rank_filter = set()


def _ictv_species_name_and_aliases(tid: int, fallback: str) -> tuple[str, str]:
    """Walk ``tid`` up to its species-rank ancestor; return the
    canonical scientific name plus the deduplicated alias string
    (``"EBV; Epstein-Barr virus; Human herpesvirus 4; ..."``) drawn
    from both the row's own taxid and the species ancestor.
    """
    if not nodes or not sci_name_by_tid:
        return fallback, ""
    species_tid = find_species_taxid(tid, nodes)
    name = sci_name_by_tid.get(species_tid) if species_tid else None
    aliases = _format_aliases(tid, species_tid, sci_name_by_tid, alias_by_tid)
    return (name or fallback), aliases

sources_for_tid: dict[int, set[str]] = {}
names_for_tid: dict[int, str] = {}
unmapped_rows: list[tuple[int, str, str, str]] = []


def _record(tid: int, name: str, source: str) -> None:
    rank = nodes.get(tid, (0, "unknown"))[1] if nodes else "unknown"
    if rank in rank_filter:
        return
    if tid in parquet_tax_ids:
        sources_for_tid.setdefault(tid, set()).add(source)
        if name and tid not in names_for_tid:
            names_for_tid[tid] = name
        return
    if params.coverage_genus_walkup and nodes:
        genus_tid = find_genus_taxid(tid, nodes)
        if genus_tid:
            rep = genus_tid if genus_tid in parquet_tax_ids else genus_to_rep_taxid.get(
                genus_tid, 0
            )
            if rep:
                sources_for_tid.setdefault(rep, set()).add(f"{source}->genus")
                if name and rep not in names_for_tid:
                    names_for_tid[rep] = name
                return
    unmapped_rows.append((tid, name, source, "absent_from_parquet"))


if "KRAKEN" in params.coverage_sources:
    kraken_df = pd.read_csv(input_.kraken_csv)
    top_kraken = (
        kraken_df.loc[kraken_df.domain == "Viruses"]
        .sort_values("percent", ascending=False)
        .head(int(params.coverage_top_n))
    )
    for r in top_kraken.itertuples():
        try:
            tid = int(r.taxonomy_id)
        except (ValueError, TypeError):
            continue
        _record(tid, str(getattr(r, "name", "")), "kraken")

if "KAIJU" in params.coverage_sources:
    try:
        kaiju_df = pd.read_csv(input_.kaiju_table, sep="\t")
    except Exception:  # noqa: BLE001
        kaiju_df = pd.DataFrame()
    if not kaiju_df.empty and "taxon_id" in kaiju_df.columns:
        kaiju_df = kaiju_df.dropna(subset=["taxon_id"])
        kaiju_df = kaiju_df.loc[kaiju_df["taxon_name"].fillna("") != "unclassified"]
        if "percent" in kaiju_df.columns:
            kaiju_df = kaiju_df.sort_values("percent", ascending=False)
        for r in kaiju_df.head(int(params.coverage_top_n)).itertuples():
            try:
                tid = int(r.taxon_id)
            except (ValueError, TypeError):
                continue
            _record(tid, str(getattr(r, "taxon_name", "")), "kaiju")

if "BLAST" in params.coverage_sources:
    for csv in input_.blastn_csvs:
        if not Path(csv).exists() or Path(csv).stat().st_size == 0:
            continue
        try:
            blast_df = pd.read_csv(csv)
        except Exception:  # noqa: BLE001
            continue
        if "accession" not in blast_df.columns:
            continue
        seen_per_csv: set[int] = set()
        for r in blast_df.itertuples():
            accession = getattr(r, "accession", None)
            if accession is None or pd.isna(accession):
                continue
            acc_str = str(accession).strip()
            tid = acc_to_tax.get(acc_str) or acc_to_tax.get(acc_str.split(".")[0])
            if tid is None:
                continue
            if tid in seen_per_csv:
                continue
            seen_per_csv.add(tid)
            _record(tid, str(getattr(r, "match_name", "")), "blast")

selected_viruses = virus_db_df[
    virus_db_df["tax_id"].astype(int).isin(sources_for_tid.keys())
]

Path(output.virus_fasta).parent.mkdir(parents=True, exist_ok=True)
with open(output.virus_fasta, "w") as f, open(output.virus_names, "w") as nf:
    nf.write("chrom\ttax_id\tname\taliases\tsources\n")
    for row in selected_viruses.itertuples():
        accession = row.name.strip().split()[0]
        tid = int(row.tax_id)
        species, aliases = _ictv_species_name_and_aliases(
            tid, names_for_tid.get(tid, "")
        )
        # Strip tabs from aliases so the TSV stays well-formed.
        aliases = aliases.replace("\t", " ")
        source_tag = ";".join(sorted(sources_for_tid.get(tid, set())))
        f.write(f">{row.name.strip()}\n{row.sequence}\n")
        nf.write(f"{accession}\t{tid}\t{species}\t{aliases}\t{source_tag}\n")

with open(output.unmapped_taxids, "w") as uf:
    uf.write("tax_id\tname\tsource\treason\n")
    for tid, name, source, reason in unmapped_rows:
        clean_name = name.replace("\t", " ").replace("\n", " ")
        uf.write(f"{tid}\t{clean_name}\t{source}\t{reason}\n")

if Path(output.virus_fasta).stat().st_size == 0:
    with open(output.virus_fasta, "w") as f:
        f.write(">DUMMY_REF\n")
        f.write("N" * 100 + "\n")

index_prefix = params.index_prefix
Path(index_prefix).parent.mkdir(parents=True, exist_ok=True)
# Shell-quote every interpolated path so a RESULTS_FOLDER / reference path
# with a space (common on macOS external volumes) does not split into
# multiple shell arguments. The cleanup glob keeps its trailing `*`
# outside the quotes so it still expands.
_q_index = shlex.quote(str(index_prefix))
_q_fasta = shlex.quote(str(output.virus_fasta))
_q_log = shlex.quote(str(log_path))
_q_r1 = shlex.quote(str(input_.r1))
_q_r2 = shlex.quote(str(input_.r2))
_q_bam = shlex.quote(str(output.bam))
_shell(f"bwa index -p {_q_index} {_q_fasta} > {_q_log} 2>&1")
_shell(
    f"bwa mem -t {threads} {_q_index} {_q_r1} {_q_r2} "
    f"| samtools sort -o {_q_bam} - >> {_q_log} 2>&1"
)
_shell(f"samtools index {_q_bam}")
_shell(f"rm -rf {_q_index}*")
