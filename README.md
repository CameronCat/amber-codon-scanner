# mcr-ptm-annotator

[![CI](https://github.com/CameronPiepkorn/mcr-ptm-annotator/actions/workflows/ci.yml/badge.svg)](https://github.com/CameronPiepkorn/mcr-ptm-annotator/actions)
[![Python 3.9+](https://img.shields.io/badge/python-3.9+-blue.svg)](https://www.python.org/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)

A Python library for annotating McrA protein sequences with known
post-translational modification (PTM) sites in the active site of
methyl-coenzyme M reductase.

Built to support research from the Nayak lab:

> Nayak, Metcalf et al. (2017). *Post-translational thioamidation of
> methyl-coenzyme M reductase, a key enzyme in methanogenic and
> methanotrophic Archaea.* eLife 6:e29218.

> Nayak, Mahanta et al. (2020). *Post-translational thioamidation of
> methyl-coenzyme M reductase is required for its activity in
> Methanosarcina acetivorans.*  PLoS Biology 18:e3000507.

---

## Background

The α subunit of MCR (McrA) contains six experimentally characterised
post-translational modifications. The first five were identified in the
crystal structure of *Methanothermobacter marburgensis* MCR
(PDB [1MRO](https://www.rcsb.org/structure/1MRO), Ermler et al. 1997).
A sixth (didehydroaspartate) was identified by mass spectrometry and
high-resolution X-ray crystallography (PDB [5A0Y](https://www.rcsb.org/structure/5A0Y),
Wagner et al. 2016):

| PTM name | Residue | Position (M. marburgensis) | PDB | Reference |
|---|---|---|---|---|
| Thioglycine | G | 445 | 1MRO (GL3) | Nayak 2017 eLife |
| Didehydroaspartate | D | 446 | 5A0Y | Wagner 2016 Angew Chem |
| 1-N-methylhistidine | H | 257 | 1MRO (MHS) | Ermler 1997 |
| 5-(S)-methylarginine | R | 271 | 1MRO (AGM) | Deobald 2018 Sci Rep |
| 2-(S)-methylglutamine | Q | 400 | 1MRO (MGN) | Ermler 1997 |
| S-methylcysteine | C | 452 | 1MRO (SMC) | Ermler 1997 |

**Known variability across organisms** (Kahnt 2007, Selmer 2000):
- Thioglycine and 1-N-methylhistidine: conserved in all methanogens examined
- 5-(S)-methylarginine: present in all methanogens examined, absent in ANME-1
- 2-(S)-methylglutamine: absent in *M. barkeri*
- S-methylcysteine: low abundance or absent in many methanogens including *M. maripaludis*
- Didehydroaspartate: present in *M. marburgensis* and *M. barkeri*, absent in *M. wolfeii*

> ⚠ **Important limitations:**
> - This tool reports **candidate** sites only. PTM presence must be
>   confirmed experimentally (mass spectrometry).
> - Position mapping is a **heuristic** scaled from PDB 1MRO. For
>   definitive mapping, align your sequence against PDB 1MRO chain A
>   (see below).
> - PTM patterns vary across archaeal lineages. ANME MCRs have a
>   different modification pattern (Shima et al. 2012).

---

## Installation

```bash
git clone https://github.com/CameronPiepkorn/mcr-ptm-annotator
cd mcr-ptm-annotator
pip install -e ".[dev]"
```

---

## Quick Start

```python
from mcr_ptm_annotator import McrAPTMAnnotator, report

annotator = McrAPTMAnnotator(
    position_window=30,        # ±residues to search around expected position
    require_residue_match=True # only report correct residue type
)

# Annotate a single sequence (McrA protein, single-letter amino acids)
hits = annotator.annotate_sequence(my_mcra_protein_sequence, seq_id="MA0528")
print(report.summary(hits, seq_id="MA0528"))

# Export
report.to_tsv(hits, seq_id="MA0528", path="ptm_hits.tsv")
report.to_json(hits, seq_id="MA0528", path="ptm_hits.json")
```

### Scan a whole FASTA file

```python
results = annotator.annotate_fasta("my_mcrA_sequences.faa")
for seq_id, hits in results.items():
    print(report.summary(hits, seq_id=seq_id))
```

### Inspect the PTM database

```python
from mcr_ptm_annotator import KNOWN_PTMS

for ptm in KNOWN_PTMS:
    print(ptm.name, ptm.position_marburgensis, ptm.references)
```

---

## How Position Mapping Works

Because McrA is highly conserved across methanogens, the known PTM
residues fall in narrow, predictable regions relative to total sequence
length. The tool scales each reference position linearly:

```
expected_pos = round(ref_pos × query_length / 553)
```

Then it searches ±30 residues around that expected position for the
correct residue type.

**This is a shortcut, not a substitute for alignment.** Confidence
levels reflect distance from the expected position:

| Confidence | Δ from expected |
|---|---|
| `high` | ≤ 10 residues |
| `moderate` | 11–20 residues |
| `low` | 21–30 residues |

---

## Definitive Position Mapping (Recommended)

For publication-quality results, align your McrA sequence against
PDB 1MRO chain A using MUSCLE or MAFFT:

```bash
# 1. Download PDB 1MRO chain A sequence
efetch -db protein -id 1MRO_A -format fasta > 1MRO_A.faa

# 2. Combine with your query sequence
cat 1MRO_A.faa my_mcrA.faa > combined.faa

# 3. Align
muscle -in combined.faa -out aligned.faa
# or
mafft --auto combined.faa > aligned.faa

# 4. Read off the alignment columns for positions 257, 271, 400, 445, 452
```

---

## Obtaining Real McrA Sequences

The example FASTA contains a synthetic placeholder. Download real
sequences from:

| Organism | Protein | Accession |
|---|---|---|
| *M. acetivorans* C2A | McrA | [NP_618892.1](https://www.ncbi.nlm.nih.gov/protein/NP_618892.1) |
| *M. marburgensis* (reference) | McrA | [UniProt P11558](https://www.uniprot.org/uniprot/P11558) |
| *M. mazei* | McrA | [NP_632996.1](https://www.ncbi.nlm.nih.gov/protein/NP_632996.1) |

```bash
efetch -db protein -id NP_618892.1 -format fasta > McrA_acetivorans.faa
```

---

## Repository Structure

```
mcr-ptm-annotator/
├── mcr_ptm_annotator/
│   ├── __init__.py
│   ├── annotator.py      ← McrAPTMAnnotator, PTMHit
│   ├── ptm_database.py   ← KNOWN_PTMS (experimentally verified only)
│   ├── report.py         ← TSV / JSON / text export
│   └── utils.py          ← FASTA parser
├── tests/
│   └── test_annotator.py
├── examples/
│   ├── example_mcrA.fasta    ← synthetic placeholder + NCBI links
│   └── basic_usage.py
├── .github/workflows/ci.yml
├── pyproject.toml
└── README.md
```

---

## Running Tests

```bash
pytest
pytest --cov=mcr_ptm_annotator
```

---

## References

- Ermler et al. (1997) *Crystal structure of methyl-coenzyme M reductase.*
  Science 278:1457. PDB: 1MRO.
- Selmer et al. (2000) *Biosynthesis of methylated amino acids in the active
  site region of MCR.* J Biol Chem 275:3755.
- Kahnt et al. (2007) *Post-translational modifications in the active site
  region of MCR from methanogenic and methanotrophic archaea.*
  FEBS J 274:4913.
- Nayak et al. (2017) *Post-translational thioamidation of MCR.*
  eLife 6:e29218.
- Wagner et al. (2016) *Didehydroaspartate modification in MCR.*
  Angew Chem Int Ed Engl 55:10630. PDB: 5A0Y.
- Deobald et al. (2018) *Radical SAM methyltransferase for
  sp3-C-methylation of arginine in MCR.* Sci Rep 8:7404.
- Nayak et al. (2020) *Thioamidation required for MCR activity.*
  PLoS Biol 18:e3000507.

---

## Citation

If you use this tool in published research, please cite:

```bibtex
@article{nayak2017,
  title   = {Post-translational thioamidation of methyl-coenzyme M
             reductase, a key enzyme in methanogenic and methanotrophic Archaea},
  author  = {Nayak, Dipti D and Mahanta, Nilkamal and Mitchell, Douglas A
             and Metcalf, William W},
  journal = {eLife},
  volume  = {6},
  pages   = {e29218},
  year    = {2017}
}
```

---

## License

MIT © 2024. See [LICENSE](LICENSE).
