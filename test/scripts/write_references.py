"""
Emit FASTA files for the synthetic host and viral references shared by the
fixture builders. Lets shell tools (bwa index, makeblastdb, kraken2-build,
kaiju-mkbwt) consume the same sequences the FASTQ synthesizer used.
"""
from __future__ import annotations

import argparse
from pathlib import Path

from synthesize_fastq import get_host_reference, get_virus_reference


def fasta_block(name: str, seq: str, line_width: int = 70) -> str:
    lines = [f">{name}"]
    for i in range(0, len(seq), line_width):
        lines.append(seq[i : i + line_width])
    return "\n".join(lines) + "\n"


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

    args.host_fasta.parent.mkdir(parents=True, exist_ok=True)
    args.host_fasta.write_text(fasta_block("synthetic_host", get_host_reference()))

    args.virus_fasta.parent.mkdir(parents=True, exist_ok=True)
    args.virus_fasta.write_text(fasta_block("synthetic_virus", get_virus_reference()))

    if args.virus_protein_fasta is not None:
        # Crude six-frame-ish protein for Kaiju. Take the nucleotide sequence
        # and translate frame 0; if any stop codons appear, just emit ASCII
        # padding so kaiju-mkbwt has something to ingest.
        nt = get_virus_reference()
        table = {
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
        protein = "".join(
            table.get(nt[i : i + 3], "X") for i in range(0, len(nt) - 2, 3)
        ).replace("*", "X")
        args.virus_protein_fasta.parent.mkdir(parents=True, exist_ok=True)
        args.virus_protein_fasta.write_text(
            fasta_block("synthetic_virus_aa", protein)
        )


if __name__ == "__main__":
    main()
