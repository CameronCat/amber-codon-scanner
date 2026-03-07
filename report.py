"""
mcr_ptm_annotator.report
~~~~~~~~~~~~~~~~~~~~~~~~
Export PTMHit results to TSV, JSON, or human-readable summary.
"""

from __future__ import annotations

import json
from pathlib import Path
from typing import Iterable

from .annotator import PTMHit
from .ptm_database import KNOWN_PTMS


TSV_HEADER = "\t".join([
    "seq_id", "ptm_name", "residue_type", "query_position",
    "residue_found", "expected_position", "position_delta",
    "confidence", "modification", "references",
])


def to_tsv(hits: Iterable[PTMHit], seq_id: str = "input",
           path: str | Path | None = None) -> str:
    rows = [TSV_HEADER]
    for h in hits:
        rows.append("\t".join([
            seq_id,
            h.ptm.name,
            h.ptm.residue_type,
            str(h.query_position),
            h.residue,
            str(h.expected_position),
            str(h.position_delta),
            h.confidence,
            h.ptm.modification,
            "; ".join(h.ptm.references),
        ]))
    content = "\n".join(rows) + "\n"
    if path:
        Path(path).write_text(content)
        return str(path)
    return content


def to_json(hits: Iterable[PTMHit], seq_id: str = "input",
            path: str | Path | None = None) -> str:
    records = [
        {
            "seq_id": seq_id,
            "ptm_name": h.ptm.name,
            "residue_type": h.ptm.residue_type,
            "query_position": h.query_position,
            "residue_found": h.residue,
            "expected_position": h.expected_position,
            "position_delta": h.position_delta,
            "confidence": h.confidence,
            "modification": h.ptm.modification,
            "pdb_code_1MRO": h.ptm.pdb_code,
            "references": list(h.ptm.references),
            "notes": h.notes,
        }
        for h in hits
    ]
    content = json.dumps(records, indent=2)
    if path:
        Path(path).write_text(content)
        return str(path)
    return content


def summary(hits: list[PTMHit], seq_id: str = "input") -> str:
    lines = [
        "=" * 74,
        f"  MCR PTM Annotator — Results for: {seq_id}",
        "=" * 74,
        f"  PTM sites found : {len(hits)} / {len(KNOWN_PTMS)} known sites",
        "=" * 74,
    ]

    if hits:
        lines.append(
            f"{'PTM name':<25} {'pos':>6} {'aa':<4} {'Δpos':>5} "
            f"{'conf':<10} PDB code"
        )
        lines.append("-" * 74)
        for h in hits:
            lines.append(
                f"{h.ptm.name:<25} {h.query_position:>6} {h.residue:<4} "
                f"{h.position_delta:>+5} {h.confidence:<10} {h.ptm.pdb_code}"
            )
    else:
        lines.append("  No PTM sites identified.")
        lines.append("  Is this sequence McrA? Expected length ~480-560 aa.")

    lines.append("=" * 74)
    lines.append(
        "\n⚠ Positions are heuristic estimates scaled from PDB 1MRO "
        "(M. marburgensis).\n"
        "  For definitive mapping, align against PDB 1MRO chain A "
        "using MUSCLE or MAFFT.\n"
        "  PTM presence requires experimental confirmation (mass spectrometry)."
    )
    return "\n".join(lines)
