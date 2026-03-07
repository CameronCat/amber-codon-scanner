"""
amber_codon_scanner.scanner
~~~~~~~~~~~~~~~~~~~~~~~~~~~
Scan archaeal DNA sequences for UAG (amber) codons and classify each
as a likely pyrrolysine-coding codon or a true translation stop.

Biological background
---------------------
Pyrrolysine (Pyl, single-letter O) is the 22nd genetically encoded amino
acid, incorporated co-translationally at UAG codons in certain methanogenic
archaea, including Methanosarcina acetivorans, M. mazei, M. barkeri, and
Methanomassiliicoccus luminyensis.

Two evidence types are used for classification:

1. PYLIS element (Pyrrolysine Insertion Sequence)
   A GC-rich stem-loop structure located downstream of UAG codons in
   pyrrolysine-encoding genes that promotes read-through by Pyl-tRNA.
   This tool searches for a GC-rich stem-loop pattern in the 80 nt
   downstream of each UAG. The detection parameters (window size,
   GC threshold, stem length) are heuristic choices — they are NOT
   calibrated from a published dataset.  Results are labelled
   "PYLIS-like" only, not "PYLIS confirmed".

2. ORF context
   A UAG that falls mid-ORF (i.e., sufficient in-frame coding sequence
   upstream, and a downstream in-frame stop within a search window) is
   a stronger pyrrolysine candidate than a terminal UAG.
   This is a necessary-but-not-sufficient heuristic with no quantitative
   calibration.

What this tool does NOT do
--------------------------
- It does not predict pyrrolysine with validated sensitivity/specificity.
  No such classifier has been published for all archaeal lineages.
- It does not perform homology search (BLAST/HMM) against pylS/pylB/pylTSBCD.
  For definitive classification, check for the pyl gene cluster in the
  same genome (see README for instructions using hmmer/NCBI).
- It does not handle frameshifts or introns.

References
----------
Shalvarjian, Chadwick et al. (2025) PNAS 122:e2517473122
  -- Nayak lab: ambiguous amber codon usage in pyrrolysine-encoding methanogens
  -- This is the primary paper this tool is designed to support.
"""

from __future__ import annotations

import re
from dataclasses import dataclass, field
from typing import Iterator

from .utils import parse_fasta, reverse_complement, translate_codon


# ── PYLIS stem-loop detection ──────────────────────────────────────────────
# The PYLIS element is a GC-rich stem-loop found downstream of UAG codons
# in pyrrolysine-encoding genes.  The detection parameters below are
# heuristic choices, NOT derived from a calibrated dataset:
#   - 80 nt search window
#   - GC-rich stem ≥5 nt + loop 4-8 nt + GC-rich run ≥4 nt
#   - Minimum 55% GC in the window
# These are intentionally permissive; results are flagged as "PYLIS-like"
# only.  Do not interpret a positive hit as confirmed pyrrolysine encoding.
_PYLIS_WINDOW = 80       # nt downstream to search
_PYLIS_STEM_MIN = 5      # minimum stem length (bp)

# GC-rich stem ≥5 nt + loop 4-8 nt + GC-rich run ≥4 nt (relaxed match)
_PYLIS_RE = re.compile(
    r"([GC]{5,})"          # stem half 1
    r"[ACGT]{4,8}"         # loop
    r"([GC]{4,})",         # stem half 2 (reverse complement not enforced,
                           # but GC-richness is a necessary condition)
    re.IGNORECASE,
)

# Minimum GC fraction in the PYLIS window to report a hit
_PYLIS_WINDOW_GC_MIN = 0.55


# ── Classification labels ──────────────────────────────────────────────────
class Classification:
    PYL_CANDIDATE   = "pyrrolysine_candidate"
    STOP_LIKELY     = "stop_likely"
    AMBIGUOUS       = "ambiguous"


@dataclass
class AmberCodon:
    """A single UAG codon found in an input sequence."""

    seq_id: str
    position: int          # 0-based nt position of the U in UAG
    strand: str            # '+' or '-'
    frame: int             # reading frame 0, 1, or 2
    mid_orf: bool          # True if UAG is flanked by in-frame codons
                           # with a downstream in-frame stop
    pylis_detected: bool   # True if PYLIS-like element found downstream
    downstream_gc: float   # GC fraction of the _PYLIS_WINDOW nt downstream
    classification: str    # Classification label (see Classification)
    notes: list[str] = field(default_factory=list)

    def __repr__(self) -> str:
        return (
            f"<AmberCodon {self.seq_id} pos={self.position}{self.strand} "
            f"frame={self.frame} class={self.classification}>"
        )


# ── Main scanner ───────────────────────────────────────────────────────────
class AmberCodonScanner:
    """
    Scan DNA sequences for UAG codons and classify them.

    Parameters
    ----------
    min_orf_length : int
        Minimum number of codons upstream of UAG (within the same frame)
        required to consider a UAG as mid-ORF (default 10).
    downstream_stop_window : int
        How many codons downstream to search for an in-frame stop codon
        that would confirm the UAG is internal to an ORF (default 50).
    check_pylis : bool
        Search for PYLIS-like stem-loop elements downstream of each UAG
        (default True).
    """

    def __init__(
        self,
        min_orf_length: int = 10,
        downstream_stop_window: int = 50,
        check_pylis: bool = True,
    ) -> None:
        self.min_orf_length = min_orf_length
        self.downstream_stop_window = downstream_stop_window
        self.check_pylis = check_pylis

    # ── Public API ──────────────────────────────────────────────────────

    def scan_sequence(
        self,
        sequence: str,
        seq_id: str = "input",
    ) -> list[AmberCodon]:
        """
        Find and classify all UAG codons in *sequence* (both strands).

        Parameters
        ----------
        sequence : str
            Raw DNA sequence (IUPAC, case-insensitive).
        seq_id : str
            Label used in AmberCodon.seq_id.

        Returns
        -------
        list[AmberCodon] sorted by position.
        """
        sequence = sequence.upper().replace(" ", "").replace("\n", "")
        results: list[AmberCodon] = []

        for codon in self._find_uag_codons(sequence, seq_id, strand="+"):
            results.append(codon)

        rc = reverse_complement(sequence)
        for codon in self._find_uag_codons(rc, seq_id, strand="-"):
            # Map position back to forward-strand coordinates
            codon.position = len(sequence) - codon.position - 3
            results.append(codon)

        results.sort(key=lambda c: c.position)
        return results

    def scan_fasta(
        self,
        fasta_path: str,
    ) -> dict[str, list[AmberCodon]]:
        """
        Scan every sequence in a FASTA file.

        Returns
        -------
        dict mapping seq_id → list[AmberCodon]
        """
        results: dict[str, list[AmberCodon]] = {}
        for seq_id, sequence in parse_fasta(fasta_path):
            results[seq_id] = self.scan_sequence(sequence, seq_id=seq_id)
        return results

    # ── Private helpers ─────────────────────────────────────────────────

    def _find_uag_codons(
        self, sequence: str, seq_id: str, strand: str
    ) -> Iterator[AmberCodon]:
        """Yield classified AmberCodon objects for all UAGs on one strand."""
        for match in re.finditer(r"(?=(TAG))", sequence, re.IGNORECASE):
            pos = match.start()
            frame = pos % 3

            mid_orf = self._is_mid_orf(sequence, pos, frame)
            pylis = False
            ds_gc = 0.0

            if self.check_pylis:
                pylis, ds_gc = self._check_pylis(sequence, pos)

            classification, notes = self._classify(mid_orf, pylis, ds_gc)

            yield AmberCodon(
                seq_id=seq_id,
                position=pos,
                strand=strand,
                frame=frame,
                mid_orf=mid_orf,
                pylis_detected=pylis,
                downstream_gc=ds_gc,
                classification=classification,
                notes=notes,
            )

    def _is_mid_orf(self, sequence: str, uag_pos: int, frame: int) -> bool:
        """
        Return True if the UAG appears to be internal to an ORF.

        Criteria (both must be met):
        1. At least min_orf_length in-frame codons upstream without a stop.
        2. At least one in-frame non-stop codon downstream before the next
           in-frame stop (within downstream_stop_window codons).
        """
        # --- Upstream: check for start and absence of stops ---
        upstream_codons = 0
        pos = uag_pos - 3
        while pos >= frame and upstream_codons < self.min_orf_length:
            codon = sequence[pos: pos + 3]
            if len(codon) < 3:
                break
            aa = translate_codon(codon)
            if aa == "*":
                return False   # hit a stop before reaching min_orf_length
            upstream_codons += 1
            pos -= 3
        if upstream_codons < self.min_orf_length:
            return False

        # --- Downstream: look for an in-frame stop within window ---
        found_downstream_stop = False
        pos = uag_pos + 3   # first codon after UAG
        codons_checked = 0
        while pos + 3 <= len(sequence) and codons_checked < self.downstream_stop_window:
            codon = sequence[pos: pos + 3]
            aa = translate_codon(codon)
            if aa == "*":
                found_downstream_stop = True
                break
            codons_checked += 1
            pos += 3

        return found_downstream_stop

    def _check_pylis(
        self, sequence: str, uag_pos: int
    ) -> tuple[bool, float]:
        """
        Search for a PYLIS-like element in the _PYLIS_WINDOW nt downstream
        of the UAG.

        Returns (pylis_detected: bool, downstream_gc: float).
        """
        start = uag_pos + 3
        end = min(len(sequence), start + _PYLIS_WINDOW)
        window = sequence[start:end]

        if not window:
            return False, 0.0

        gc = (window.count("G") + window.count("C")) / len(window)

        if gc < _PYLIS_WINDOW_GC_MIN:
            return False, gc

        # Check for stem-loop pattern
        if _PYLIS_RE.search(window):
            return True, gc

        return False, gc

    def _classify(
        self, mid_orf: bool, pylis: bool, ds_gc: float
    ) -> tuple[str, list[str]]:
        """
        Apply classification rules.

        Rules (conservative — designed to minimise false positives):
        - PYLIS detected AND mid-ORF  → pyrrolysine_candidate (strong)
        - PYLIS detected, not mid-ORF → ambiguous
        - mid-ORF only, no PYLIS       → ambiguous
        - neither                      → stop_likely
        """
        notes: list[str] = []

        if pylis and mid_orf:
            notes.append("PYLIS-like element detected downstream")
            notes.append("UAG is internal to a predicted ORF")
            return Classification.PYL_CANDIDATE, notes

        if pylis and not mid_orf:
            notes.append("PYLIS-like element detected but UAG is not clearly mid-ORF")
            return Classification.AMBIGUOUS, notes

        if mid_orf and not pylis:
            notes.append("UAG is internal to a predicted ORF but no PYLIS element found")
            notes.append("Consider checking for pylTSBCD gene cluster in genome")
            return Classification.AMBIGUOUS, notes

        notes.append("No evidence for pyrrolysine recoding")
        return Classification.STOP_LIKELY, notes
