# functions.py

from pathlib import Path
import re
import subprocess
import os
import pandas as pd
import numpy as np
import pyfastx
import janitor  # pandas-flavor extension
import seaborn as sns
import matplotlib.pyplot as plt
import panel as pn
import altair as alt
from bs4 import BeautifulSoup
import warnings
from io import StringIO

warnings.filterwarnings("ignore")


def read_file_as_blob(file_path: str) -> str:
    """
    Reads a file and returns its content in hexadecimal format.

    Args:
        file_path (str): Path to the file to read.

    Returns:
        str: Hexadecimal representation of the file content.
    """
    with open(file_path, 'rb') as file:
        blob_data = file.read().hex()
    return blob_data


def common_suffix(folder: str) -> str:
    """
    Determines the common suffix among sequencing files in a given folder.

    Args:
        folder (str): Path to the folder containing sequencing files.

    Returns:
        str: The common suffix shared by the sequencing files.
    """
    samples = sorted(
        [file.name for file in Path(folder).iterdir() if re.search(r"\.(fq|fastq|fa|fasta|fna)$", file.name)]
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
    """
    Identifies paired-end read prefixes from sequencing files in a folder.

    Args:
        folder (str): Path to the folder containing sequencing files.

    Returns:
        list: A list of common prefixes for paired-end reads.
    """
    def common_prefix(str1: str, str2: str) -> str:
        prefix = ""
        for a, b in zip(str1, str2):
            if a != b:
                break
            prefix += a
        return prefix

    samples = sorted(
        [x.stem for x in Path(folder).iterdir() if re.search(r"\.(fq|fastq|fa|fasta|fna)$", x.name)]
    )

    prefixes = []
    for i in range(0, len(samples), 2):
        read1, read2 = samples[i], samples[i + 1]
        common = common_prefix(read1, read2)
        prefixes.append(common)

    return prefixes


def kaiju_db_files(kaiju_db: str) -> tuple:
    """
    Retrieves necessary Kaiju database files from a given directory.

    Args:
        kaiju_db (str): Path to the Kaiju database directory.

    Returns:
        tuple: Paths to the .fmi, names.dmp, and nodes.dmp files.
    """
    files = [x for x in Path(kaiju_db).iterdir() if x.is_file()]
    fmi = next((x for x in files if x.suffix == ".fmi"), None)
    names = next((x for x in files if x.name == "names.dmp"), None)
    nodes = next((x for x in files if x.name == "nodes.dmp"), None)
    return fmi, names, nodes


def fastx_file_to_df(fastx_file: str) -> pd.DataFrame:
    """
    Converts a FASTA/FASTQ file into a pandas DataFrame.

    Args:
        fastx_file (str): Path to the FASTA/FASTQ file.

    Returns:
        pd.DataFrame: DataFrame containing sequence names and sequences.
    """
    fastx = pyfastx.Fastx(fastx_file)
    reads = list(zip(*[(entry.name, entry.sequence) for entry in fastx]))

    df = (
        pd.DataFrame({"name": reads[0], "sequence": reads[1]})
        .assign(read_len=lambda x: x.sequence.str.len())
        .sort_values("read_len", ascending=False)
    )

    return df


def wrangle_kraken(kraken_file: str) -> pd.DataFrame:
    """
    Processes Kraken2 output into a structured pandas DataFrame.

    Args:
        kraken_file (str): Path to the Kraken2 output file.

    Returns:
        pd.DataFrame: Processed Kraken2 data.
    """
    kraken = (
        pd.read_csv(
            kraken_file, sep="\t", header=None,
            names=["percent", "count_clades", "count", "tax_lvl", "taxonomy_id", "name"]
        )
        .assign(name=lambda x: x.name.str.strip())
        .assign(
            domain=lambda x: np.select(
                [x.tax_lvl.isin(["D", "U", "R"])],
                [x.name],
                default=pd.NA
            )
        )
        .fillna(method="ffill")
    )

    return kraken


def run_blastn(contigs_csv: str, db: str, temp_file: str, threads: int) -> pd.DataFrame:
    """
    Runs BLASTN on contigs to find matches in a given database.

    Args:
        contigs_csv (str): Path to the CSV file containing contig sequences.
        db (str): Path to the BLASTN database.
        temp_file (str): Temporary file to store individual contig sequences.
        threads (int): Number of threads to use for BLASTN.

    Returns:
        pd.DataFrame: DataFrame with BLASTN results appended to the contig data.
    """
    # Ensure BLASTDB environment variable is set
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
            "-outfmt", "6 stitle sacc pident slen"
        ]
        match = subprocess.check_output(command, universal_newlines=True).strip()
        matches.append(match)

    df = df.assign(matches=matches).loc[lambda x: x.matches != ""]
    if df.empty:
        return df

    df[["match_name", "accession", "percent_identity", "sequence_len"]] = (
        df.matches.str.split("\t", expand=True).iloc[:, :4]
    )
    df = df.assign(sequence_len=lambda x: x.sequence_len.str.split("\n").str[0])

    return df


def parse_bwa_flagstat(flagstat_file: str) -> tuple:
    """
    Parses BWA flagstat output to extract total reads and percentage mapped.

    Args:
        flagstat_file (str): Path to the BWA flagstat output file.

    Returns:
        tuple: Total number of reads and percentage of reads mapped.
    """
    pattern_total = r"(\d+) \+ \d+ paired in sequencing"
    pattern_mapped = r"(\d+) \+ \d+ with itself and mate mapped"

    with open(flagstat_file) as f:
        flagstat = f.read()

    total_reads = int(re.search(pattern_total, flagstat).group(1))
    total_mapped = int(re.search(pattern_mapped, flagstat).group(1))

    percent_mapped = (total_mapped / total_reads) * 100 if total_reads > 0 else 0

    return total_reads, percent_mapped


def parse_fastp(fastp_report: str) -> pd.DataFrame:
    """
    Extracts summary statistics from a FASTP HTML report.

    Args:
        fastp_report (str): Path to the FASTP HTML report.

    Returns:
        pd.DataFrame: DataFrame containing summary statistics.
    """
    summary_info = {}

    with open(fastp_report, "r") as f:
        soup = BeautifulSoup(f, "html.parser")

    for table in soup.find_all("table", class_="summary_table")[:4]:
        for row in table.find_all("tr"):
            cells = row.find_all("td")
            if len(cells) >= 2:
                key = cells[0].text.strip()
                value = cells[1].text.strip()
                summary_info[key] = value

    df = pd.DataFrame.from_dict(summary_info, orient="index", columns=["value"])
    df = df.rename_axis("description").reset_index()

    return df


def plot_flagstat(flagstat_file: str) -> alt.Chart:
    """
    Generates an Altair plot showing the proportion of reads aligned vs. unaligned.

    Args:
        flagstat_file (str): Path to the BWA flagstat output file.

    Returns:
        alt.Chart: Altair chart object.
    """
    total_reads, percent_aligned = parse_bwa_flagstat(flagstat_file)
    number_aligned = int(total_reads * percent_aligned / 100)
    number_unaligned = total_reads - number_aligned

    # Create DataFrame
    df = pd.DataFrame(
        {"amount": [number_unaligned, number_aligned], "type": ["Unaligned", "Aligned"]}
    )

    # Generate plot
    plot = (
        alt.Chart(df)
        .mark_bar()
        .encode(
            x=alt.X(
                "sum(amount)",
                stack="normalize",
                axis=alt.Axis(format="%"),
                title=None,
            ),
            color=alt.Color("type:N", scale=alt.Scale(scheme='dark2')),
            tooltip=[
                alt.Tooltip("amount:Q", title="Number of reads"),
                alt.Tooltip("type:N"),
            ],
        )
        .properties(
            title="Reads Aligned to Host",
            width="container",
            height=100,
        )
    )

    return plot


def plot_kaiju(kaiju_table: str, cutoff: float = 0.01, max_entries: int = 10) -> alt.Chart:
    """
    Generates an Altair plot for Kaiju classification results.

    Args:
        kaiju_table (str): Path to the Kaiju output table.
        cutoff (float, optional): Minimum percentage to include a taxon. Defaults to 0.01.
        max_entries (int, optional): Maximum number of taxa to display. Defaults to 10.

    Returns:
        alt.Chart: Altair chart object.
    """
    df = pd.read_csv(kaiju_table, sep="\t")
    df = df.assign(percent=lambda x: x.percent / 100)

    unclassified_percent = df.loc[df.taxon_name == "unclassified", "percent"].squeeze()
    df = (
        df[df.taxon_name != "unclassified"]
        .query("percent > @cutoff")
        .nlargest(max_entries, "percent")
    )

    plot = (
        alt.Chart(df)
        .mark_bar()
        .encode(
            x=alt.X(
                "percent:Q",
                axis=alt.Axis(format="%"),
                title=f"Percent of Reads ({unclassified_percent*100:.1f}% Unclassified)",
                scale=alt.Scale(zero=True)
            ),
            y=alt.Y("taxon_name:N", sort="-x", title=None),
            color=alt.Color("taxon_name:N", legend=None),
            tooltip=[
                alt.Tooltip("taxon_name:N", title="Taxon"),
                alt.Tooltip("reads:Q", title="Number of Reads")
            ],
        )
        .properties(
            width="container",
            title="Kaiju Classification"
        )
    )

    return plot


def kraken_df(
    kraken_csv: str,
    level: str = "species",
    cutoff: float = 0.01,
    max_entries: int = 10,
    virus_only: bool = True,
) -> tuple:
    """
    Processes Kraken2 classification results into a DataFrame.

    Args:
        kraken_csv (str): Path to the processed Kraken2 CSV file.
        level (str, optional): Taxonomic level to filter. Defaults to "species".
        cutoff (float, optional): Minimum percentage to include a taxon. Defaults to 0.01.
        max_entries (int, optional): Maximum number of taxa to display. Defaults to 10.
        virus_only (bool, optional): Whether to include only viral taxa. Defaults to True.

    Returns:
        tuple: DataFrame of filtered taxa and the percentage unclassified.
    """
    taxonomy = {
        "domain": "D",
        "phylum": "P",
        "class": "C",
        "order": "O",
        "family": "F",
        "genus": "G",
        "species": "S",
    }

    df = pd.read_csv(kraken_csv)
    df = df.assign(percent=lambda x: x.percent / 100)

    unclassified_percent = df[df.domain == "unclassified"]["percent"].sum()

    df = df[df.tax_lvl == taxonomy.get(level, "S")]
    df = df[df.percent > cutoff]
    df = df.sort_values("percent", ascending=False)

    if virus_only:
        df = df[df.domain == "Viruses"]

    df = df.head(max_entries)

    return df, unclassified_percent


def plot_kraken(
    kraken_csv: str,
    level: str = "species",
    cutoff: float = 0.001,
    max_entries: int = 10,
    virus_only: bool = True,
) -> alt.Chart:
    """
    Generates an Altair plot for Kraken2 classification results.

    Args:
        kraken_csv (str): Path to the processed Kraken2 CSV file.
        level (str, optional): Taxonomic level to filter. Defaults to "species".
        cutoff (float, optional): Minimum percentage to include a taxon. Defaults to 0.001.
        max_entries (int, optional): Maximum number of taxa to display. Defaults to 10.
        virus_only (bool, optional): Whether to include only viral taxa. Defaults to True.

    Returns:
        alt.Chart: Altair chart object.
    """
    df, unclassified_percent = kraken_df(
        kraken_csv, level, cutoff, max_entries, virus_only
    )

    plot = (
        alt.Chart(df)
        .mark_bar()
        .encode(
            x=alt.X(
                "percent:Q",
                axis=alt.Axis(format="%"),
                title=f"Percent of Reads ({unclassified_percent*100:.1f}% Unclassified)"
            ),
            y=alt.Y("name:N", sort="-x", title=None),
            color=alt.Color("name:N", legend=None),
            tooltip=[
                alt.Tooltip("domain:N", title="Domain"),
                alt.Tooltip("count_clades:Q", title="Number of Reads")
            ],
        )
        .properties(
            width="container",
            title="Kraken Classification"
        )
    )

    return plot


def plot_blastn(blastn_csv: str) -> alt.Chart:
    """
    Generates an Altair plot for BLASTN annotation results of contigs.

    Args:
        blastn_csv (str): Path to the BLASTN results CSV file.

    Returns:
        alt.Chart: Altair chart object.
    """
    df = pd.read_csv(blastn_csv)

    if df.empty:
        df = pd.DataFrame({"match_name": ["No Classified Contigs"], "count": [0]})

    df_grouped = df.groupby("match_name").size().reset_index(name='count')

    plot = (
        alt.Chart(df_grouped)
        .mark_bar()
        .encode(
            x=alt.X(
                "count:Q",
                title="Number of Contigs",
                axis=alt.Axis(tickMinStep=1),
            ),
            y=alt.Y("match_name:N", sort="-x", title=None),
            color=alt.Color("match_name:N", legend=None),
            tooltip=[
                alt.Tooltip("match_name:N", title="Match Name"),
                alt.Tooltip("count:Q", title="Number of Contigs")
            ],
        )
        .properties(
            width="container",
            title="BLASTN Annotation of Contigs"
        )
    )

    return plot


def alignment_stats(flagstat_file: str, species: str) -> tuple:
    """
    Compiles alignment statistics into markdown and plots.

    Args:
        flagstat_file (str): Path to the BWA flagstat output file.
        species (str): Name of the species the reads were aligned to.

    Returns:
        tuple: Markdown pane and Vega pane for alignment stats.
    """
    total_reads, percent_aligned = parse_bwa_flagstat(flagstat_file)
    number_aligned = int(total_reads * percent_aligned / 100)
    number_unaligned = total_reads - number_aligned

    # Markdown summary
    alignment_summary = pn.pane.Markdown(
        f"""
        ### Total Number of Reads:
        {total_reads:,}
        ### Reads Aligned to {species} Genome:
        {number_aligned:,} ({percent_aligned:.2f}%)
        ### Reads NOT Aligned to {species} Genome:
        {number_unaligned:,} ({100 - percent_aligned:.2f}%)
        """,
        name=f"{species} Alignment Stats"
    )

    # Alignment plot
    flagstat_plot = plot_flagstat(flagstat_file).interactive()
    flagstat_pane = pn.pane.Vega(
        flagstat_plot, sizing_mode="stretch_both", name=f"{species} Alignment Plot"
    )

    return alignment_summary, flagstat_pane


def panel_report(
    result_folder: str,
    blastn_file: str,
    kraken_file: str,
    kaiju_table: str,
    secondary_host: str = None
) -> pn.Column:
    """
    Generates an interactive HTML report using Panel.

    Args:
        result_folder (str): Path to the result folder for the sample.
        blastn_file (str): Path to the BLASTN results CSV file.
        kraken_file (str): Path to the processed Kraken2 CSV file.
        kaiju_table (str): Path to the Kaiju output table.
        secondary_host (str, optional): Name of the secondary host species. Defaults to None.

    Returns:
        pn.Column: Panel object representing the report.
    """
    result_folder = Path(result_folder)
    sample_name = result_folder.name

    # Initialize Panel extensions
    pn.extension("tabulator")
    pn.extension("vega", sizing_mode="stretch_width", template="fast")
    pn.widgets.Tabulator.theme = 'modern'

    def header(
        text: str,
        bg_color: str = "#04c273",
        height: int = 150,
        fontsize: str = "20px",
        textalign: str = "center"
    ) -> pn.pane.Markdown:
        """
        Creates a styled header for the report sections.

        Args:
            text (str): Header text.
            bg_color (str, optional): Background color. Defaults to "#04c273".
            height (int, optional): Height of the header. Defaults to 150.
            fontsize (str, optional): Font size. Defaults to "20px".
            textalign (str, optional): Text alignment. Defaults to "center".

        Returns:
            pn.pane.Markdown: Styled markdown pane.
        """
        return pn.pane.Markdown(
            text,
            styles={
                "color": "white",
                "padding": "10px",
                "text-align": textalign,
                "font-size": fontsize,
                "background": bg_color,
                "margin": "10px",
                "height": f"{height}px",
            }
        )

    # Alignment and Read Statistics
    flagstat_human_log = result_folder / "logs" / "human_contamination_flagstat.txt"
    human_stats, human_pane = alignment_stats(str(flagstat_human_log), species="Human")

    # Optional secondary host alignment
    if secondary_host:
        flagstat_secondary_log = result_folder / "logs" / "secondary_contamination_flagstat.txt"
        if flagstat_secondary_log.exists():
            secondary_stats, secondary_pane = alignment_stats(
                str(flagstat_secondary_log), species=secondary_host
            )
        else:
            secondary_stats, secondary_pane = None, None
    else:
        secondary_stats, secondary_pane = None, None

    # FASTP report
    fastp_report = next(result_folder.rglob("FASTP/*.html"), None)
    if fastp_report:
        fastp_df = parse_fastp(str(fastp_report))
        fastp_table = pn.widgets.Tabulator(
            fastp_df,
            layout='fit_columns',
            show_index=False,
            name="Read Summary from FASTP"
        )
    else:
        fastp_table = pn.pane.Markdown("No FASTP report found.")

    # Alignment section header
    alignment_subheader = header(
        text="""
        ## Alignment and Read Statistics
        Reads were aligned to Human and optionally other host species using BWA.
        """,
        bg_color="#04c273",
        height=80,
        textalign="left"
    )

    # Alignment stats tabs
    alignment_tabs = [("Human Alignment", pn.Column(human_stats, pn.layout.Divider(), human_pane))]
    if secondary_stats:
        alignment_tabs.append(
            (f"{secondary_host} Alignment", pn.Column(secondary_stats, pn.layout.Divider(), secondary_pane))
        )
    alignment_tabs.append(("FASTP Summary", fastp_table))

    alignment_section = pn.Column(alignment_subheader, pn.Tabs(*alignment_tabs))

    # Classification of Raw Reads
    kraken_plot = plot_kraken(kraken_file).interactive()
    kraken_domain_plot = plot_kraken(kraken_file, level="domain", virus_only=False).interactive()
    kaiju_plot = plot_kaiju(kaiju_table).interactive()

    kraken_pane = pn.pane.Vega(kraken_plot, sizing_mode="stretch_both", name="Kraken (Viruses)")
    kraken_domain_pane = pn.pane.Vega(kraken_domain_plot, sizing_mode="stretch_both", name="Kraken (All Domains)")
    kaiju_pane = pn.pane.Vega(kaiju_plot, sizing_mode="stretch_both", name="Kaiju")

    raw_header = header(
        text="""
        ## Classification of Raw Reads
        Reads were classified using Kaiju and Kraken2.
        """,
        bg_color="#04c273",
        height=80,
        textalign="left"
    )

    raw_section = pn.Column(raw_header, pn.Tabs(kraken_pane, kraken_domain_pane, kaiju_pane))

    # Contig Classification
    blastn_plot = plot_blastn(blastn_file).interactive()
    blastn_pane = pn.pane.Vega(blastn_plot, sizing_mode="stretch_both", name="BLASTN")

    blastn_df = pd.read_csv(blastn_file)
    if blastn_df.empty:
        blastn_df = pd.DataFrame({"sequence": ["NO SEQUENCES GENERATED"]})
    else:
        blastn_df = blastn_df.drop(columns=["name", "matches"])

    blastn_table = pn.widgets.Tabulator(
        blastn_df,
        layout='fit_columns',
        pagination='local',
        page_size=15,
        show_index=False,
        name="Contig Table"
    )

    contig_header = header(
        text="""
        ## Classification of Contigs
        Contigs assembled using MEGAHIT and annotated with BLASTN.
        """,
        bg_color="#04c273",
        height=120,
        textalign="left"
    )

    contig_section = pn.Column(contig_header, pn.Tabs(blastn_pane, blastn_table))

    # Coverage Plots
    coverage_plot_path = result_folder / "COVERAGE_PLOTS"
    coverage_plots = list(coverage_plot_path.glob("*.svg"))

    coverage_tabs = []
    if coverage_plots:
        for plot_file in coverage_plots:
            name = plot_file.stem[:20]
            svg_pane = pn.pane.SVG(plot_file, sizing_mode='stretch_width', name=name)
            coverage_tabs.append((name, svg_pane))
    else:
        coverage_tabs.append(("No Coverage Plots", pn.pane.Markdown("## No Coverage Plots Available")))

    coverage_header = header(
        text="## Alignment Coverage",
        bg_color="#04c273",
        height=80,
        textalign="left"
    )

    coverage_section = pn.Column(coverage_header, pn.Tabs(*coverage_tabs))

    # Assemble the report
    main_header = header(
        text=f"""
        # Virushanter Report
        ## Sample: {sample_name}
        """,
        fontsize="24px",
        bg_color="#011a01",
        height=185
    )

    all_tabs = pn.Tabs(
        ("Alignment Stats", alignment_section),
        ("Classification of Raw Reads", raw_section),
        ("Classification of Contigs", contig_section),
        ("Alignment Coverage", coverage_section),
        tabs_location="left",
    )

    report = pn.Column(
        main_header,
        pn.layout.Divider(),
        all_tabs,
    )

    return report
