# provenance.smk
#
# Capture what actually ran: the conda-resolved version of every tool the
# pipeline invoked, plus a machine-readable run provenance sidecar the
# HTML report renders.
#
# A probe rule that references the same envs/<env>.yaml as the real work
# reuses the identical materialised conda prefix under --sdm conda, so its
# conda-meta directory lists the exact resolved name-version-build. The
# probe dumps only the filename stems (no Python needed, since a minimal
# tool env such as fastp carries no interpreter); the stems are parsed in
# scripts/collect_software_versions.py.

import platform as _platform

_BATCH = Path(config["SAMPLES"]).name

# Envs a default run always exercises.
_VERSION_ENVS = [
    "fastp",
    "bwa",
    "samtools",
    "kraken",
    "kaiju",
    "megahit",
    "spades",
    "pilon",
    "blastn",
    "checkv",
    "mosdepth",
]
# Optional envs, probed only when their stage runs so software_versions
# records exactly what executed (no version for a tool that never ran).
if config.get("MULTIQC", "TRUE") == "TRUE":
    _VERSION_ENVS.append("multiqc")
if config.get("GENOMAD", "FALSE") == "TRUE":
    _VERSION_ENVS.append("genomad")
if config.get("QUAST", "FALSE") == "TRUE":
    _VERSION_ENVS.append("quast")
if config.get("HOST_REMOVAL", "bwa") == "hostile":
    _VERSION_ENVS.append("hostile")


def _reporthanter_pin() -> str:
    """Return the reportHanter git ref pinned in envs/reporthanter.yaml.

    The report renderer is pip-installed from a tag, so it is absent from
    conda-meta and the probe cannot see it; the configured tag is the
    honest provenance value. Returns an empty string if unparsable.
    """
    env_path = Path(workflow.basedir) / "envs" / "reporthanter.yaml"
    try:
        for line in env_path.read_text().splitlines():
            marker = "reporthanter.git@"
            if marker in line:
                return line.split(marker, 1)[1].strip()
    except OSError:
        pass
    return ""


def _snakemake_version() -> str:
    try:
        import snakemake

        return getattr(snakemake, "__version__", "")
    except Exception:  # noqa: BLE001
        return ""


# Rule: probe one conda env for its resolved package versions.
rule software_version_probe:
    output:
        tsv=f"{RESULT_FOLDER}/logs/versions/{{env}}.tsv",
    params:
        env=lambda wildcards: wildcards.env,
    wildcard_constraints:
        env="[A-Za-z0-9]+",
    conda:
        lambda wildcards: f"../envs/{wildcards.env}.yaml"
    resources:
        mem_mb=500,
        runtime=5,
    shell:
        r"""
        mkdir -p "$(dirname {output.tsv})"
        : > {output.tsv}
        for f in "$CONDA_PREFIX"/conda-meta/*.json; do
            [ -e "$f" ] || continue
            b=${{f##*/}}
            printf '%s\t%s\n' "{params.env}" "${{b%.json}}" >> {output.tsv}
        done
        """


# Rule: merge every per-env probe into one resolved-version table.
rule collect_software_versions:
    input:
        expand(f"{RESULT_FOLDER}/logs/versions/{{env}}.tsv", env=_VERSION_ENVS),
    output:
        tsv=f"{RESULT_FOLDER}/software_versions.tsv",
    resources:
        mem_mb=1000,
        runtime=10,
    conda:
        "../envs/panel.yaml"
    script:
        "../scripts/collect_software_versions.py"


# Rule: write the per-run provenance sidecar (DB build identity + tool
# versions) that the report renders and operators diff across runs.
rule write_provenance:
    input:
        software_versions=rules.collect_software_versions.output.tsv,
    output:
        json=f"{RESULT_FOLDER}/run_provenance_{_BATCH}.json",
        tsv=f"{RESULT_FOLDER}/run_provenance_{_BATCH}.tsv",
    params:
        run_name=_BATCH,
        assemblers=ASSEMBLERS,
        reporthanter_version=_reporthanter_pin(),
        snakemake_version=_snakemake_version(),
        python_version=_platform.python_version(),
    resources:
        mem_mb=1000,
        runtime=10,
    conda:
        "../envs/panel.yaml"
    script:
        "../scripts/write_provenance.py"
