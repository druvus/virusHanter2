"""Snakemake script: aggregate Kaiju output to a per-taxon table,
then rewrite ``taxon_name`` with the ICTV-binomial species name.

Runs ``kaiju2table`` in the Kaiju conda env (so the binary is on
PATH), then applies the same species-rank walk-up the BLAST
canonicaliser uses so the Classification of reads tab and the
Dashboard Kaiju card show the same canonical species names as
the Assembly classification chart.
"""

import subprocess
import sys
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from scripts.functions import canonicalise_taxon_names  # noqa: E402

snakemake = snakemake  # type: ignore[name-defined]

input_ = snakemake.input
output = snakemake.output
params = snakemake.params

subprocess.run(
    [
        "kaiju2table",
        "-t", str(input_.nodes),
        "-n", str(input_.names),
        "-r", "genus",
        "-e",
        "-o", str(output.kaiju_table),
        str(input_.kaiju_out),
    ],
    check=True,
)

# Canonicalise taxon_name -> ICTV-binomial species via the published
# taxdump alongside the parquet. Falls back to a no-op when
# TAXDUMP_NODES is unset / missing - the kaiju2table output is then
# left as kaiju2table produced it.
taxdump_nodes = str(params.taxdump_nodes) if params.taxdump_nodes else ""
names_dmp = ""
if taxdump_nodes:
    names_dmp_path = Path(taxdump_nodes).parent / "names.dmp"
    if names_dmp_path.is_file():
        names_dmp = str(names_dmp_path)

if taxdump_nodes and names_dmp:
    df = pd.read_csv(output.kaiju_table, sep="\t")
    df = canonicalise_taxon_names(
        df,
        taxid_col="taxon_id",
        name_col="taxon_name",
        nodes_dmp=taxdump_nodes,
        names_dmp=names_dmp,
    )
    df.to_csv(output.kaiju_table, sep="\t", index=False)
