
from __future__ import annotations

from pathlib import Path
from typing import Iterator


_COMPLEMENT = str.maketrans("ACGTacgt", "TGCAtgca")


def reverse_complement(seq: str) -> str:
    return seq.translate(_COMPLEMENT)[::-1]


def parse_fasta(path: str | Path) -> Iterator[tuple[str, str]]:
    path = Path(path)
    current_id: str | None = None
    buffer: list[str] = []

    with path.open() as fh:
        for raw_line in fh:
            line = raw_line.strip()
            if not line or line.startswith(";"):
                continue
            if line.startswith(">"):
                if current_id is not None:
                    yield current_id, "".join(buffer)
                current_id = line[1:].split()[0]
                buffer = []
            else:
                buffer.append(line.upper())

    if current_id is not None:
        yield current_id, "".join(buffer)


_CODON_TABLE: dict[str, str] = {
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


def translate_codon(codon: str) -> str:
    return _CODON_TABLE.get(codon.upper(), "X")
