"""
mcr_ptm_annotator.annotator
~~~~~~~~~~~~~~~~~~~~~~~~~~~
Annotate McrA protein sequences with known post-translational modification
sites, mapped by position range scanning.

⚠ IMPORTANT — WHAT THIS TOOL DOES AND DOES NOT DO
---------------------------------------------------
This tool identifies candidate residues in a query McrA sequence that
correspond to known PTM sites in the M. marburgensis reference (PDB 1MRO).

It does NOT:
  - Confirm that a modification is actually present (mass spec is required)
  - Perform full sequence alignment (use MUSCLE/MAFFT for that)
  - Predict novel PTMs not listed in the database

Mapping approach
----------------
Because McrA is highly conserved across methanogens, the known PTM residues
(G, H, R, Q, C) are found in narrow, predictable position windows relative
to the total sequence length. We scan a ±window around the scaled expected
position and report candidate residues of the correct amino acid type.

This is a HEURISTIC shortcut. For definitive mapping, align your sequence
against PDB 1MRO chain A using MUSCLE or MAFFT and read off the alignment
column. Instructions for this are in the README.

References
----------
Ermler et al. (1997) Science 278:1457  -- PDB 1MRO structure, 5 PTMs
Wagner et al. (2016) Angew Chem 55:10630 -- didehydroaspartate, 6th PTM (PDB 5A0Y)
Nayak et al. (2017) eLife 6:e29218     -- thioamidation genetics
Nayak et al. (2020) PLoS Biol 18:e3000507 -- PTM functional interactions
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Optional

from .ptm_database import KNOWN_PTMS, KnownPTM
from .utils import parse_fasta


# Reference sequence length for M. marburgensis McrA (UniProt P11558)
# Used to scale expected positions to query sequence length.
_MARBURGENSIS_MCRA_LENGTH = 553

# ±window (residues) to search around the scaled expected position.
# Chosen to be wide enough to accommodate typical sequence length variation
# across methanogens (~10% length difference) while remaining specific.
_POSITION_WINDOW = 30  # residues


@dataclass
class PTMHit:
    """A candidate PTM site identified in a query sequence."""

    ptm: KnownPTM
    query_position: int        # 1-based position in query sequence
    residue: str               # amino acid at this position (should match ptm.residue_type)
    expected_position: int     # scaled expected position in query
    position_delta: int        # query_position - expected_position
    confidence: str            # 'high' | 'moderate' | 'low'
    notes: list[str] = field(default_factory=list)

    def __repr__(self) -> str:
        return (
            f"<PTMHit {self.ptm.name} at pos={self.query_position} "
            f"({self.residue}) conf={self.confidence}>"
        )


class McrAPTMAnnotator:
    """
    Annotate an McrA protein sequence with candidate PTM sites.

    Parameters
    ----------
    position_window : int
        ±window in residues around the scaled expected position to search
        (default 30). Increase for highly divergent sequences.
    require_residue_match : bool
        If True (default), only report hits where the residue type matches
        the known PTM residue type. Set to False for exploratory analysis.
    """

    def __init__(
        self,
        position_window: int = _POSITION_WINDOW,
        require_residue_match: bool = True,
    ) -> None:
        self.position_window = position_window
        self.require_residue_match = require_residue_match

    # ── Public API ────────────────────────────────────────────────────────

    def annotate_sequence(
        self,
        sequence: str,
        seq_id: str = "input",
    ) -> list[PTMHit]:
        """
        Find candidate PTM sites in a single McrA protein sequence.

        Parameters
        ----------
        sequence : str
            Single-letter amino acid sequence (case-insensitive).
            Should be the McrA (alpha subunit) sequence only.
        seq_id : str
            Label for reporting.

        Returns
        -------
        list[PTMHit] sorted by query position.
        """
        sequence = sequence.upper().strip()
        seq_len = len(sequence)
        hits: list[PTMHit] = []

        for ptm in KNOWN_PTMS:
            expected = self._scale_position(
                ptm.position_marburgensis, seq_len
            )
            candidates = self._find_candidates(sequence, ptm, expected)
            hits.extend(candidates)

        hits.sort(key=lambda h: h.query_position)
        return hits

    def annotate_fasta(
        self,
        fasta_path: str,
    ) -> dict[str, list[PTMHit]]:
        """
        Annotate every sequence in a FASTA file.

        Returns
        -------
        dict mapping seq_id → list[PTMHit]
        """
        results: dict[str, list[PTMHit]] = {}
        for seq_id, sequence in parse_fasta(fasta_path):
            results[seq_id] = self.annotate_sequence(sequence, seq_id=seq_id)
        return results

    # ── Private helpers ───────────────────────────────────────────────────

    def _scale_position(self, ref_position: int, query_length: int) -> int:
        """
        Scale a M. marburgensis reference position to the query sequence length.

        Uses simple linear scaling:
            expected = round(ref_position * query_length / ref_length)
        """
        scaled = round(ref_position * query_length / _MARBURGENSIS_MCRA_LENGTH)
        return max(1, min(scaled, query_length))

    def _find_candidates(
        self,
        sequence: str,
        ptm: KnownPTM,
        expected_pos: int,
    ) -> list[PTMHit]:
        """
        Search ±window around expected_pos for residues matching ptm.residue_type.
        Returns a list of PTMHit objects (usually 0 or 1).
        """
        lo = max(0, expected_pos - self.position_window - 1)
        hi = min(len(sequence), expected_pos + self.position_window)
        window_seq = sequence[lo:hi]

        hits: list[PTMHit] = []
        for i, aa in enumerate(window_seq):
            query_pos = lo + i + 1  # 1-based
            if self.require_residue_match and aa != ptm.residue_type:
                continue

            delta = query_pos - expected_pos
            confidence = self._confidence(delta)

            notes = [
                f"Reference position: {ptm.residue_type}{ptm.position_marburgensis} "
                f"in M. marburgensis (PDB 1MRO)",
                f"Scaled expected position in query: {expected_pos}",
            ]
            if abs(delta) > 15:
                notes.append(
                    "Position is >15 residues from expected — verify by "
                    "alignment against PDB 1MRO chain A"
                )

            hits.append(PTMHit(
                ptm=ptm,
                query_position=query_pos,
                residue=aa,
                expected_position=expected_pos,
                position_delta=delta,
                confidence=confidence,
                notes=notes,
            ))

        # If multiple candidates, return only the closest to expected position
        if len(hits) > 1:
            hits = [min(hits, key=lambda h: abs(h.position_delta))]

        return hits

    def _confidence(self, delta: int) -> str:
        """
        Assign confidence based on distance from expected position.

        'high'     : within ±10 residues (typical for closely related MCRs)
        'moderate' : within ±20 residues
        'low'      : within ±window but >20 residues away
        """
        abs_delta = abs(delta)
        if abs_delta <= 10:
            return "high"
        elif abs_delta <= 20:
            return "moderate"
        else:
            return "low"
