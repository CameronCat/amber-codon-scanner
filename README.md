# amber-codon-scanner

[![CI](https://github.com/CameronPiepkorn/amber-codon-scanner/actions/workflows/ci.yml/badge.svg)](https://github.com/CameronPiepkorn/amber-codon-scanner/actions)
[![Python 3.9+](https://img.shields.io/badge/python-3.9+-blue.svg)](https://www.python.org/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)

A Python library for scanning archaeal DNA sequences for UAG (amber) codons
and classifying each as a likely **pyrrolysine-coding** codon or a true
**translation stop**.

Built to support research from the Nayak lab:

> Shalvarjian, Chadwick et al. (2025). *Methanogenic archaea encoding
> pyrrolysine maintain ambiguous amber codon usage.*
> PNAS 122(45):e2517473122.

---

## Background

Pyrrolysine (Pyl) is the 22nd genetically encoded amino acid, incorporated
at UAG codons in certain methanogenic archaea including *Methanosarcina
acetivorans*, *M. mazei*, and *M. barkeri*. These organisms use UAG
ambiguously — the same codon serves as both a stop signal and a
pyrrolysine codon depending on context.

This tool uses two evidence sources to classify each UAG:

| Evidence | Description | Source |
|---|---|---|
| **Mid-ORF context** | UAG flanked by in-frame coding sequence with a downstream in-frame stop | Heuristic |
| **PYLIS element** | GC-rich stem-loop downstream that promotes Pyl-tRNA read-through | Heuristic pattern (see PYLIS section below) |

> ⚠ **Important:** Classification is heuristic only. No validated
> sensitivity/specificity data exist for this classifier across all archaeal
> lineages. For definitive results, check for the **pylTSBCD gene cluster**
> in the genome (see below) and/or use experimental data.

---

## Installation

```bash
# From source
git clone https://github.com/CameronPiepkorn/amber-codon-scanner
cd amber-codon-scanner
pip install -e ".[dev]"
```

---

## Quick Start

```python
from amber_codon_scanner import AmberCodonScanner, Classification, report

scanner = AmberCodonScanner(
    min_orf_length=10,          # codons upstream required to call mid-ORF
    downstream_stop_window=50,  # codons downstream to search for in-frame stop
    check_pylis=True,           # search for PYLIS-like element
)

codons = scanner.scan_sequence(my_dna_string, seq_id="MA0859_mttB")
print(report.summary(codons))

# Export
report.to_tsv(codons, "results.tsv")
report.to_json(codons, "results.json")

# Filter to candidates only
candidates = [c for c in codons if c.classification == Classification.PYL_CANDIDATE]
```

### Scan a whole FASTA file

```python
results = scanner.scan_fasta("my_genome.fasta")
for seq_id, codons in results.items():
    print(f"\n=== {seq_id} ===")
    print(report.summary(codons))
```

---

## Classification Rules

| mid-ORF | PYLIS detected | Classification |
|---|---|---|
| ✅ | ✅ | `pyrrolysine_candidate` |
| ✅ | ❌ | `ambiguous` |
| ❌ | ✅ | `ambiguous` |
| ❌ | ❌ | `stop_likely` |

---

## Confirming Results — pylTSBCD Gene Cluster

The most reliable confirmation of pyrrolysine encoding is the presence of
the **pyl biosynthesis gene cluster** in the same genome:

| Gene | Function |
|---|---|
| `pylS` | Pyrrolysyl-tRNA synthetase (charges tRNA with pyrrolysine) |
| `pylT` | Pyrrolysine tRNA (anticodon CUA, decodes UAG) |
| `pylB/C/D` | Pyrrolysine biosynthesis enzymes |

To check for these genes using HMMER against your genome:

```bash
# Download pyl HMM profiles from Pfam / TIGRFAM
# Then:
hmmsearch --tblout pyl_hits.txt Pyl_synthetase.hmm my_genome.faa
```

Known pyrrolysine-encoding genes in *M. acetivorans* C2A (GenBank AE010299.1):

| Locus | Gene | Contains UAG |
|---|---|---|
| MA0859 | mttB | Yes (trimethylamine methyltransferase) |
| MA4384 | mtmB | Yes (monomethylamine methyltransferase) |

For the full list of pyrrolysine-encoding methyltransferases and their
locus tags, query GenBank AE010299.1 directly or see the NCBI gene pages
linked in the example FASTA file header.

---

## PYLIS Element Detection

The PYLIS (Pyrrolysine Insertion Sequence) element is a GC-rich stem-loop
found downstream of UAG codons in pyrrolysine-encoding genes. It promotes
read-through of UAG by the pyrrolysyl-tRNA.

This tool uses a relaxed heuristic pattern (GC-rich stem ≥5 bp + loop
4-8 nt + GC-rich continuation, within an 80 nt window downstream of UAG)
to flag PYLIS-like structures. **The detection parameters are heuristic
choices, not calibrated from a published benchmark.** This is intentionally
permissive and will produce false positives in GC-rich genomes.
Results are labelled `pylis_detected: True` only, not "PYLIS confirmed".

---

## Repository Structure

```
amber-codon-scanner/
├── amber_codon_scanner/
│   ├── __init__.py
│   ├── scanner.py      ← AmberCodonScanner, AmberCodon, Classification
│   ├── report.py       ← TSV / JSON / text export
│   └── utils.py        ← FASTA parser, codon table, reverse complement
├── tests/
│   └── test_scanner.py
├── examples/
│   ├── example_sequences.fasta   ← synthetic placeholders + NCBI links
│   └── basic_usage.py
├── .github/workflows/ci.yml
├── pyproject.toml
└── README.md
```

---

## Running Tests

```bash
pytest
pytest --cov=amber_codon_scanner
```

---

## Citation

If you use this tool in published research, please cite:

```bibtex
@article{shalvarjian2025,
  title   = {Methanogenic archaea encoding pyrrolysine maintain
             ambiguous amber codon usage},
  author  = {Shalvarjian, Katherine E and Chadwick, Garrett L and
             P{\'e}rez, Pilar I and Woods, Patrick H and Orphan, Victoria J
             and Nayak, Dipti D},
  journal = {Proceedings of the National Academy of Sciences},
  volume  = {122},
  number  = {45},
  pages   = {e2517473122},
  year    = {2025}
}
```

---

## License

MIT © 2024. See [LICENSE](LICENSE).
