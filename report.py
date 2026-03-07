"""
amber_codon_scanner.report
~~~~~~~~~~~~~~~~~~~~~~~~~~
Export AmberCodon results to TSV, JSON, or a human-readable summary.
"""

from __future__ import annotations

import json
from pathlib import Path
from typing import Iterable

from .scanner import AmberCodon, Classification


TSV_HEADER = "\t".join([
    "seq_id", "position", "strand", "frame",
    "mid_orf", "pylis_detected", "downstream_gc",
    "classification", "notes",
])


def to_tsv(codons: Iterable[AmberCodon], path: str | Path | None = None) -> str:
    rows = [TSV_HEADER]
    for c in codons:
        rows.append("\t".join([
            c.seq_id,
            str(c.position),
            c.strand,
            str(c.frame),
            str(c.mid_orf),
            str(c.pylis_detected),
            f"{c.downstream_gc:.3f}",
            c.classification,
            "; ".join(c.notes),
        ]))
    content = "\n".join(rows) + "\n"
    if path:
        Path(path).write_text(content)
        return str(path)
    return content


def to_json(codons: Iterable[AmberCodon], path: str | Path | None = None) -> str:
    records = [
        {
            "seq_id": c.seq_id,
            "position": c.position,
            "strand": c.strand,
            "frame": c.frame,
            "mid_orf": c.mid_orf,
            "pylis_detected": c.pylis_detected,
            "downstream_gc": round(c.downstream_gc, 4),
            "classification": c.classification,
            "notes": c.notes,
        }
        for c in codons
    ]
    content = json.dumps(records, indent=2)
    if path:
        Path(path).write_text(content)
        return str(path)
    return content


def summary(codons: list[AmberCodon], show_all: bool = False) -> str:
    """
    Return a formatted text summary.

    If show_all is False (default), only pyrrolysine_candidate and
    ambiguous codons are shown.  Set show_all=True to include stop_likely.
    """
    if not show_all:
        display = [c for c in codons if c.classification != Classification.STOP_LIKELY]
    else:
        display = list(codons)

    total = len(codons)
    n_pyl  = sum(1 for c in codons if c.classification == Classification.PYL_CANDIDATE)
    n_amb  = sum(1 for c in codons if c.classification == Classification.AMBIGUOUS)
    n_stop = sum(1 for c in codons if c.classification == Classification.STOP_LIKELY)

    lines = [
        "=" * 74,
        "  Amber Codon Scanner — Results",
        "=" * 74,
        f"  Total UAG codons found : {total}",
        f"  Pyrrolysine candidates : {n_pyl}",
        f"  Ambiguous              : {n_amb}",
        f"  Stop (likely)          : {n_stop}",
        "=" * 74,
    ]

    if display:
        lines.append(
            f"{'seq_id':<20} {'pos':>8} {'str':<3} {'fr':<3} "
            f"{'mid_orf':<8} {'PYLIS':<6} {'class':<25} notes"
        )
        lines.append("-" * 74)
        for c in display:
            note_str = c.notes[0] if c.notes else ""
            lines.append(
                f"{c.seq_id[:20]:<20} {c.position:>8} {c.strand:<3} {c.frame:<3} "
                f"{str(c.mid_orf):<8} {str(c.pylis_detected):<6} "
                f"{c.classification:<25} {note_str}"
            )

    lines.append("=" * 74)
    lines.append(
        "\n⚠ Classification is heuristic only. Confirm with pylTSBCD gene\n"
        "  cluster presence and/or experimental data. See README for details."
    )
    return "\n".join(lines)
