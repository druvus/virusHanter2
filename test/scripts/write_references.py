"""
Emit FASTA files for the synthetic host and viral references shared by the
fixture builders. Lets shell tools (bwa index, makeblastdb, kraken2-build,
kaiju-mkbwt) consume the same sequences the FASTQ synthesizer used.

There is one host record and three viral records (alpha / beta / gamma).
The viral nucleotide FASTA carries Kraken-style taxid tags so
``kraken2-build --add-to-library`` can pick them up; the protein FASTA is
formatted for ``kaiju-mkbwt`` (header ``>P0000N_<taxid>``).
"""
from __future__ import annotations

import argparse
from pathlib import Path

from synthesize_fastq import get_host_reference, get_virus_references


CODON_TABLE = {
    "TTT": "F", "TTC": "F", "TTA": "L", "TTG": "L",
    "CTT": "L", "CTC": "L", "CTA": "L", "CTG": "L",
    "ATT": "I", "ATC": "I", "ATA": "I", "ATG": "M",
    "GTT": "V", "GTC": "V", "GTA": "V", "GTG": "V",
    "TCT": "S", "TCC": "S", "TCA": "S", "TCG": "S",
    "CCT": "P", "CCC": "P", "CCA": "P", "CCG": "P",
    "ACT": "T", "ACC": "T", "ACA": "T", "ACG": "T",
    "GCT": "A", "GCC": "A", "GCA": "A", "GCG": "A",
    "TAT": "Y", "TAC": "Y", "TAA": "*", "TAG": "*",
    "CAT": "H", "CAC": "H", "CAA": "Q", "CAG": "Q",
    "AAT": "N", "AAC": "N", "AAA": "K", "AAG": "K",
    "GAT": "D", "GAC": "D", "GAA": "E", "GAG": "E",
    "TGT": "C", "TGC": "C", "TGA": "*", "TGG": "W",
    "CGT": "R", "CGC": "R", "CGA": "R", "CGG": "R",
    "AGT": "S", "AGC": "S", "AGA": "R", "AGG": "R",
    "GGT": "G", "GGC": "G", "GGA": "G", "GGG": "G",
}


def fasta_block(name: str, seq: str, line_width: int = 70) -> str:
    lines = [f">{name}"]
    for i in range(0, len(seq), line_width):
        lines.append(seq[i : i + line_width])
    return "\n".join(lines) + "\n"


def translate_frame0(nt: str) -> str:
    """Frame-0 translation with stop codons replaced by 'X' so the result
    is a single contiguous protein string suitable for kaiju-mkbwt.
    """
    protein = "".join(
        CODON_TABLE.get(nt[i : i + 3], "X") for i in range(0, len(nt) - 2, 3)
    )
    return protein.replace("*", "X")


def main() -> None:
    p = argparse.ArgumentParser()
    p.add_argument("--host-fasta", type=Path, required=True)
    p.add_argument("--virus-fasta", type=Path, required=True)
    p.add_argument(
        "--virus-protein-fasta",
        type=Path,
        required=False,
        help="Optional protein-translation FASTA for Kaiju mkbwt.",
    )
    args = p.parse_args()

    # Host FASTA.
    args.host_fasta.parent.mkdir(parents=True, exist_ok=True)
    args.host_fasta.write_text(fasta_block("synthetic_host", get_host_reference()))

    # Multi-record viral nucleotide FASTA with clean headers. The
    # kraken2-build step tags these with `|kraken:taxid|<id>` in its
    # own copy (see rule `kraken_library` in build_fixtures.smk); we
    # keep the canonical FASTA clean so BLAST and the virus parquet
    # see record names like "alpha" that match Kraken2's species
    # names ("synthetic virus alpha") via substring lookup downstream
    # in `attribute_contigs`.
    viruses = get_virus_references()
    args.virus_fasta.parent.mkdir(parents=True, exist_ok=True)
    with args.virus_fasta.open("w") as fh:
        for name, info in viruses.items():
            fh.write(fasta_block(name, info["sequence"]))

    if args.virus_protein_fasta is not None:
        # Kaiju ingests one protein record per virus with a sequential
        # `>P0000N_<taxid>` header. Frame-0 translation; stops -> X so
        # the protein is one contiguous string.
        args.virus_protein_fasta.parent.mkdir(parents=True, exist_ok=True)
        with args.virus_protein_fasta.open("w") as fh:
            for idx, (name, info) in enumerate(viruses.items(), start=1):
                header = f"P{idx:05d}_{info['taxid']}"
                fh.write(fasta_block(header, translate_frame0(info["sequence"])))


if __name__ == "__main__":
    main()
