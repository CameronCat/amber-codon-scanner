"""
amber_codon_scanner.utils
~~~~~~~~~~~~~~~~~~~~~~~~~
Sequence utilities shared across the package.
"""

from __future__ import annotations

from pathlib import Path
from typing import Iterator


_COMPLEMENT = str.maketrans("ACGTacgt", "TGCAtgca")


def reverse_complement(seq: str) -> str:
    """Return the reverse complement of a DNA string."""
    return seq.translate(_COMPLEMENT)[::-1]


def parse_fasta(path: str | Path) -> Iterator[tuple[str, str]]:
    """
    Minimal FASTA parser.  Yields (seq_id, sequence) tuples.
    Does not require BioPython.
    """
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


# Standard genetic code (NCBI transl_table=1).
# UAG → "*" (stop) under the standard code.
# In pyrrolysine-encoding archaea, UAG is contextually reassigned to
# pyrrolysine (single-letter O) when pylT tRNA and PYLIS element are present.
# This is NOT a formal NCBI transl_table — it is a context-dependent
# suppression, not a genome-wide code change.
# We use the standard code here and flag UAG positions separately.
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
    """
    Translate a single codon to single-letter amino acid.
    Returns 'X' for unknown/ambiguous codons.
    """
    return _CODON_TABLE.get(codon.upper(), "X")
