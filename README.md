
# Virushanter2: A Snakemake Pipeline for Viral Metagenomics Analysis

![Virushanter Logo](resources/logo.png)

## Overview

**Virushanter2** is a scalable and reproducible bioinformatics pipeline designed for the analysis of viral metagenomic sequencing data. Built using [Snakemake](https://snakemake.readthedocs.io/en/stable/), the pipeline performs quality control, host contamination removal, taxonomic classification, assembly, annotation, and report generation.

## Features

- **Quality Control**: Trimming and filtering of raw reads using `fastp`.
- **Host Contamination Removal**: Alignment to host genomes (e.g.,
  human) using `BWA` and removal of host-derived reads.
- **PCR Duplicate Accounting**: `samtools markdup -s` on the
  host-aligned BAM, useful for hybrid-capture libraries
  (Twist VRP).
- **Taxonomic Classification**: Classification of reads using `Kaiju`
  and `Kraken2`.
- **Assembly and Annotation**: Assembly of reads into contigs using
  `MEGAHIT`, polishing with `Pilon`, and annotation using `BLASTN`
  and `CheckV`.
- **Per-Reference Coverage Stats**: `mosdepth` summaries with
  thresholds at 1x/5x/10x for every reference in the Kraken-top-N
  alignment.
- **Visualization and Reporting**: Interactive per-sample HTML
  reports (Panel), per-reference SVG coverage profiles (`bam2plot`),
  and a workflow-level `MultiQC` dashboard.
- **Per-(sample, virus) detail CSV**: One row per detected Kraken
  viral taxid per sample with read counts, RPM, contig count,
  completeness at 5x, and mean coverage. Documented in
  `docs/PER_VIRUS_OUTPUT.md`.
- **Multi-run aggregation**: `scripts/merge_runs.py` combines several
  Illumina runs into master CSVs without re-running the workflow.
- **Modularity**: Organized into separate rule files for better
  maintainability.
- **Reproducibility**: Each rule specifies its own Conda environment,
  ensuring consistent computational environments across runs.

## Installation

### Prerequisites

- **Operating System**: Linux or macOS (Windows is not officially supported).
- **Conda**: [Miniconda](https://docs.conda.io/en/latest/miniconda.html) or [Anaconda](https://www.anaconda.com/products/individual) installed.
- **Snakemake**: Install via Conda:

  ```bash
  conda create -n snakemake_env -c bioconda -c conda-forge snakemake
  conda activate snakemake_env
  ```

### Clone the Repository

```bash
git clone https://github.com/yourusername/virushanter.git
cd virushanter
```


## Dependencies

All software dependencies are managed via Conda environments specified in the `envs/` directory. Snakemake will automatically create these environments when running the pipeline with the `--use-conda` flag.

## Configuration

### Setting Up `config.yaml`

Copy the example configuration file and edit it according to your setup:

```bash
cp config/config.yaml.example config/config.yaml
```

Open `config/config.yaml` and modify the following parameters:

```yaml
# Path to the directory containing your raw sequencing data
SAMPLES: "/path/to/your/samples"

# Path to the directory where results will be stored
RESULTS_FOLDER: "/path/to/your/results"

# Number of threads to use for multi-threaded tools
THREADS: 8

# Whether to clean up intermediate files after the pipeline finishes
CLEAN: "FALSE"

# Path to the BWA index for the human genome
HUMAN_INDEX: "/path/to/human/genome/index"

# Path to the Kaiju database directory
KAIJU_DB: "/path/to/kaiju/database"

# Path to the Kraken2 database
KRAKEN_DB: "/path/to/kraken2/database"

# Path to the BLASTN nucleotide database
BLASTN_DB: "/path/to/blastn/database"

# Minimum contig length to keep after assembly
CONTIG_LENGTH: 1000

# Path to the CheckV database
CHECKV_DB: "/path/to/checkv/database"

# Path to the Parquet file containing viral sequences
VIRUS_PARQUET: "/path/to/virus_sequences.parquet"

# Coverage threshold for generating coverage plots
PLOT_THRESHOLD: 10

# Number of top references to include in coverage plots
NUMBER_OF_PLOTS: 10

# Optional: Path to the BWA index for a secondary host genome
# SECONDARY_HOST_INDEX: "/path/to/secondary/host/genome/index"
# SECONDARY_HOST_NAME: "Mouse"
```

**Note**: Ensure all paths are absolute and accessible.

## Usage

### Step 1: Prepare the Environment

Activate the Snakemake environment:

```bash
conda activate snakemake_env
```

### Step 2: Configure the Pipeline

Edit the `config/config.yaml` file as described above.

### Step 3: Run the Pipeline

Execute the pipeline using the following command:

```bash
snakemake --use-conda --cores <number_of_cores>
```

Replace `<number_of_cores>` with the number of CPU cores you wish to allocate.

**Example**:

```bash
snakemake --use-conda --cores 8
```

### Dry Run

To perform a dry run and see what steps will be executed:

```bash
snakemake --use-conda --cores <number_of_cores> -n
```

### Unlock the Working Directory

If you need to unlock the working directory (e.g., after an interruption):

```bash
snakemake --unlock
```


## Pipeline Workflow

### Pre-processing

1. **Quality Control**: Trimming and filtering of raw reads using `fastp`.
2. **Host Contamination Removal**: Alignment to the human genome using `BWA` and removal of human reads.
3. **Secondary Host Removal (Optional)**: Alignment to a secondary host genome and removal of those reads.

### Classification

4. **Kaiju Classification**: Taxonomic classification of reads using `Kaiju`.
5. **Kraken2 Classification**: Taxonomic classification of reads using `Kraken2`.

### Assembly

6. **Assembly with MEGAHIT**: Assembling reads into contigs.
7. **Polishing with Pilon**: Improving assembly quality.
8. **Annotation with BLASTN**: Annotating contigs by aligning to a nucleotide database.
9. **CheckV Analysis**: Assessing viral contigs for completeness and contamination.

### Post-processing

10. **Alignment to Viral References**: Aligning reads to top viral hits for coverage analysis.
11. **Coverage Plots**: Generating coverage plots using `bam2plot`.
12. **Report Generation**: Creating interactive HTML reports using `Panel`.
13. **Run Information Aggregation**: Compiling run statistics across samples.
14. **Cleanup (Optional)**: Removing intermediate files to save space.

## Conda Environments

The `envs/` directory contains Conda environment files for each tool:

- **blastn.yaml**: For BLASTN.
- **bwa.yaml**: For BWA and SAMtools.
- **checkv.yaml**: For CheckV.
- **fastp.yaml**: For fastp.
- **kaiju.yaml**: For Kaiju.
- **kraken.yaml**: For Kraken2.
- **megahit.yaml**: For MEGAHIT.
- **panel.yaml**: For Python packages used in report generation.
- **pilon.yaml**: For Pilon.
- **samtools.yaml**: For SAMtools.

These files specify the exact versions of software to ensure reproducibility.

## Outputs

The pipeline generates the following outputs, all under
`{RESULTS_FOLDER}/{batch}/`, where `batch` is the basename of `SAMPLES`:

### Per-sample outputs (`{batch}/{sample}/`)

- **Quality Control Reports**: HTML and JSON reports from `fastp` under
  `FASTP/`.
- **Trimmed Reads**: FASTQ files after quality control (under
  `FASTP/`).
- **Host-Removed Reads**: FASTQ files with human reads removed (under
  `bwa/`).
- **PCR Duplicate Stats**: `samtools markdup -s` summary at
  `logs/human_markdup_stats.txt` (reporting only — does not filter
  reads).
- **Classification Results**:
  - Kaiju output files and classification tables under `KAIJU/`.
  - Kraken2 reports and the processed `<sample>.kraken.csv` under
    `KRAKEN/`.
- **Assembly and Annotation**:
  - Assembled contigs from MEGAHIT under `MEGAHIT/`.
  - Polished contigs from Pilon under `PILON/`.
  - BLASTN annotation results under `BLASTN/`.
  - CheckV contamination assessment under `CHECKV/`, plus the
    BLASTN+CheckV inner-joined `<sample>.merged.csv`.
- **Coverage Plots**: SVG coverage profiles for top Kraken viral hits
  under `COVERAGE_PLOTS/`.
- **Per-Reference Coverage Stats**: `mosdepth` numeric summary +
  per-region BED + thresholds (1x/5x/10x) under `MOSDEPTH/`.
- **Per-Sample HTML Report**: `REPORT/<sample>.html`.
- **Per-Sample Per-Virus CSV**: `<sample>.per_virus.csv` — one row per
  detected Kraken viral taxid with read counts, RPM, contig count,
  Kaiju match, and coverage metrics. Schema in
  `docs/PER_VIRUS_OUTPUT.md`.

### Per-batch outputs (`{batch}/`)

- **`run_information_<batch>.csv`** — one row per sample, summary
  metrics (parity-locked to the original virusHanter; see
  `docs/PARITY_NOTES.md`).
- **`per_virus_<batch>.csv`** — concatenation of every sample's
  `per_virus.csv`. The collaborator-facing detail file.
- **`multiqc_report.html`** + `multiqc_data/` — workflow-level QC
  dashboard covering fastp, samtools, kraken2, mosdepth, and markdup
  across the whole batch. Gated by `MULTIQC: "TRUE"` in
  `config.yaml` (default on).

### Combining multiple Illumina runs

`scripts/merge_runs.py` is a standalone CLI that walks several batch
folders and writes `master_per_sample.csv` and `master_per_virus.csv`
into an out-dir. Each input row already carries its `run_name`, so
the master files are self-describing.

```bash
python scripts/merge_runs.py \
    --result-folder /path/to/RESULTS/<batch1> \
    --result-folder /path/to/RESULTS/<batch2> \
    --out-dir /path/to/master/
```

See `docs/PER_VIRUS_OUTPUT.md` for the full schema and a worked
example.


## Customization

### Adding or Modifying Rules

The pipeline is modular, with rules organized in the `rules/` directory. You can add or modify rules as needed.

### Updating Conda Environments

If you need to update or add dependencies, modify the corresponding `.yaml` files in the `envs/` directory.

### Changing Parameters

Adjust parameters such as thread counts, coverage thresholds, and contig length cutoffs in the `config/config.yaml` file.

## Support

If you encounter any issues or have questions, please open an issue on the [GitHub repository](https://github.com/yourusername/virushanter/issues).

## License

This project is licensed under the MIT License. See the [LICENSE](LICENSE) file for details.
