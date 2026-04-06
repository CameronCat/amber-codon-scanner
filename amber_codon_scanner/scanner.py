
from __future__ import annotations

import re
from dataclasses import dataclass, field
from typing import Iterator

from .utils import parse_fasta, reverse_complement, translate_codon


_PYLIS_WINDOW = 80
_PYLIS_STEM_MIN = 5

_PYLIS_RE = re.compile(
    r"([GC]{5,})"
    r"[ACGT]{4,8}"
    r"([GC]{4,})",
    re.IGNORECASE,
)

_PYLIS_WINDOW_GC_MIN = 0.55


class Classification:
    PYL_CANDIDATE   = "pyrrolysine_candidate"
    STOP_LIKELY     = "stop_likely"
    AMBIGUOUS       = "ambiguous"


@dataclass
class AmberCodon:

    seq_id: str
    position: int
    strand: str
    frame: int
    mid_orf: bool
    pylis_detected: bool
    downstream_gc: float
    classification: str
    notes: list[str] = field(default_factory=list)

    def __repr__(self) -> str:
        return (
            f"<AmberCodon {self.seq_id} pos={self.position}{self.strand} "
            f"frame={self.frame} class={self.classification}>"
        )


class AmberCodonScanner:

    def __init__(
        self,
        min_orf_length: int = 10,
        downstream_stop_window: int = 50,
        check_pylis: bool = True,
    ) -> None:
        self.min_orf_length = min_orf_length
        self.downstream_stop_window = downstream_stop_window
        self.check_pylis = check_pylis


    def scan_sequence(
        self,
        sequence: str,
        seq_id: str = "input",
    ) -> list[AmberCodon]:
        sequence = sequence.upper().replace(" ", "").replace("\n", "")
        results: list[AmberCodon] = []

        for codon in self._find_uag_codons(sequence, seq_id, strand="+"):
            results.append(codon)

        rc = reverse_complement(sequence)
        for codon in self._find_uag_codons(rc, seq_id, strand="-"):
            codon.position = len(sequence) - codon.position - 3
            results.append(codon)

        results.sort(key=lambda c: c.position)
        return results

    def scan_fasta(
        self,
        fasta_path: str,
    ) -> dict[str, list[AmberCodon]]:
        results: dict[str, list[AmberCodon]] = {}
        for seq_id, sequence in parse_fasta(fasta_path):
            results[seq_id] = self.scan_sequence(sequence, seq_id=seq_id)
        return results


    def _find_uag_codons(
        self, sequence: str, seq_id: str, strand: str
    ) -> Iterator[AmberCodon]:
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
        upstream_codons = 0
        pos = uag_pos - 3
        while pos >= frame and upstream_codons < self.min_orf_length:
            codon = sequence[pos: pos + 3]
            if len(codon) < 3:
                break
            aa = translate_codon(codon)
            if aa == "*":
                return False
            upstream_codons += 1
            pos -= 3
        if upstream_codons < self.min_orf_length:
            return False

        found_downstream_stop = False
        pos = uag_pos + 3
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
        start = uag_pos + 3
        end = min(len(sequence), start + _PYLIS_WINDOW)
        window = sequence[start:end]

        if not window:
            return False, 0.0

        gc = (window.count("G") + window.count("C")) / len(window)

        if gc < _PYLIS_WINDOW_GC_MIN:
            return False, gc

        if _PYLIS_RE.search(window):
            return True, gc

        return False, gc

    def _classify(
        self, mid_orf: bool, pylis: bool, ds_gc: float
    ) -> tuple[str, list[str]]:
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
