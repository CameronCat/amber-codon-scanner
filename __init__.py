"""
amber_codon_scanner
~~~~~~~~~~~~~~~~~~~
Scan archaeal DNA sequences for UAG (amber) codons and classify each
as a likely pyrrolysine-coding codon or a true translation stop.

Based on:
  Shalvarjian, Chadwick et al. (2025) PNAS 122:e2517473122
    -- Ambiguous amber codon usage in pyrrolysine-encoding methanogens
  Namy et al. (2004) Mol Cell 13(1):1-12
    -- PYLIS element characterisation

Quick start::

    from amber_codon_scanner import AmberCodonScanner, report

    scanner = AmberCodonScanner()
    codons = scanner.scan_sequence(my_dna, seq_id="MA0528")

    print(report.summary(codons))
    report.to_tsv(codons, "results.tsv")
"""

from .scanner import AmberCodonScanner, AmberCodon, Classification
from . import report

__version__ = "0.1.0"
__all__ = ["AmberCodonScanner", "AmberCodon", "Classification", "report"]
