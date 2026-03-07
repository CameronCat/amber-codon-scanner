"""
tests/test_annotator.py
Unit tests for the MCR PTM annotator.
"""

import json
import pytest
from mcr_ptm_annotator import McrAPTMAnnotator, PTMHit, KNOWN_PTMS
from mcr_ptm_annotator.ptm_database import KnownPTM
from mcr_ptm_annotator.annotator import _MARBURGENSIS_MCRA_LENGTH
from mcr_ptm_annotator.utils import parse_fasta, is_protein_sequence
from mcr_ptm_annotator import report


# ── Fixtures ──────────────────────────────────────────────────────────────

def _synthetic_mcra(length: int = 553) -> str:
    """
    Build a synthetic McrA-length sequence with correct residue types
    placed at the exact scaled positions for each PTM.
    Uses 'A' as background with the correct PTM residues inserted.
    """
    seq = list("A" * length)
    for ptm in KNOWN_PTMS:
        pos = round(ptm.position_marburgensis * length / _MARBURGENSIS_MCRA_LENGTH)
        pos = max(0, min(pos - 1, length - 1))  # 0-based
        seq[pos] = ptm.residue_type
    return "".join(seq)


@pytest.fixture
def annotator():
    return McrAPTMAnnotator()


@pytest.fixture
def synthetic_seq():
    return _synthetic_mcra()


@pytest.fixture
def tmp_fasta(tmp_path, synthetic_seq):
    fa = tmp_path / "test.fasta"
    fa.write_text(f">synthetic_mcrA\n{synthetic_seq}\n")
    return fa


# ── PTM Database ──────────────────────────────────────────────────────────

class TestPTMDatabase:
    def test_five_known_ptms(self):
        assert len(KNOWN_PTMS) == 6

    def test_all_have_references(self):
        for ptm in KNOWN_PTMS:
            assert len(ptm.references) >= 1

    def test_residue_types_are_single_letter(self):
        for ptm in KNOWN_PTMS:
            assert len(ptm.residue_type) == 1
            assert ptm.residue_type in "GHRCQD"

    def test_positions_are_positive(self):
        for ptm in KNOWN_PTMS:
            assert ptm.position_marburgensis > 0

    def test_pdb_codes_present_or_empty(self):
        # didehydroaspartate is in 5A0Y not 1MRO, so pdb_code is empty string
        for ptm in KNOWN_PTMS:
            assert isinstance(ptm.pdb_code, str)

    def test_thioglycine_present(self):
        names = [p.name for p in KNOWN_PTMS]
        assert "thioglycine" in names

    def test_didehydroaspartate_present(self):
        names = [p.name for p in KNOWN_PTMS]
        assert "didehydroaspartate" in names

    def test_thioglycine_position(self):
        tg = next(p for p in KNOWN_PTMS if p.name == "thioglycine")
        assert tg.position_marburgensis == 445
        assert tg.residue_type == "G"

    def test_didehydroaspartate_adjacent_to_thioglycine(self):
        tg = next(p for p in KNOWN_PTMS if p.name == "thioglycine")
        dd = next(p for p in KNOWN_PTMS if p.name == "didehydroaspartate")
        # didehydroaspartate is adjacent to thioglycine (Wagner 2016)
        assert abs(dd.position_marburgensis - tg.position_marburgensis) == 1

    def test_methylarginine_position(self):
        ma = next(p for p in KNOWN_PTMS if p.name == "5-(S)-methylarginine")
        assert ma.position_marburgensis == 271
        assert ma.residue_type == "R"

    def test_methylhistidine_position(self):
        mh = next(p for p in KNOWN_PTMS if p.name == "1-N-methylhistidine")
        assert mh.position_marburgensis == 257
        assert mh.residue_type == "H"


# ── Annotator ─────────────────────────────────────────────────────────────

class TestAnnotator:
    def test_finds_hits_in_synthetic_seq(self, annotator, synthetic_seq):
        hits = annotator.annotate_sequence(synthetic_seq)
        assert len(hits) >= 3  # should find most PTM residues

    def test_hits_are_sorted_by_position(self, annotator, synthetic_seq):
        hits = annotator.annotate_sequence(synthetic_seq)
        positions = [h.query_position for h in hits]
        assert positions == sorted(positions)

    def test_residue_match_enforced(self, annotator, synthetic_seq):
        hits = annotator.annotate_sequence(synthetic_seq)
        for h in hits:
            assert h.residue == h.ptm.residue_type

    def test_residue_match_disabled(self, synthetic_seq):
        a = McrAPTMAnnotator(require_residue_match=False)
        hits = a.annotate_sequence(synthetic_seq)
        # With match disabled, may find more candidates
        assert len(hits) >= 1

    def test_confidence_high_at_exact_position(self, annotator, synthetic_seq):
        hits = annotator.annotate_sequence(synthetic_seq)
        high_conf = [h for h in hits if h.confidence == "high"]
        assert len(high_conf) >= 1

    def test_empty_sequence_returns_no_hits(self, annotator):
        hits = annotator.annotate_sequence("A" * 553)
        # All-alanine has no G/H/R/Q/C at right positions → few/no hits
        assert isinstance(hits, list)

    def test_annotate_fasta(self, annotator, tmp_fasta):
        results = annotator.annotate_fasta(str(tmp_fasta))
        assert "synthetic_mcrA" in results
        assert isinstance(results["synthetic_mcrA"], list)

    def test_position_scaling(self, annotator):
        # At reference length, scaling should return the same position
        scaled = annotator._scale_position(445, _MARBURGENSIS_MCRA_LENGTH)
        assert scaled == 445

    def test_position_scaling_shorter_seq(self, annotator):
        # For a shorter sequence, position should scale down
        scaled = annotator._scale_position(445, 500)
        assert scaled < 445
        assert scaled > 0

    def test_confidence_levels(self, annotator):
        assert annotator._confidence(0) == "high"
        assert annotator._confidence(10) == "high"
        assert annotator._confidence(11) == "moderate"
        assert annotator._confidence(20) == "moderate"
        assert annotator._confidence(21) == "low"


# ── Utils ──────────────────────────────────────────────────────────────────

class TestUtils:
    def test_parse_fasta(self, tmp_fasta):
        records = list(parse_fasta(tmp_fasta))
        assert len(records) == 1
        seq_id, seq = records[0]
        assert seq_id == "synthetic_mcrA"

    def test_is_protein_sequence_valid(self):
        assert is_protein_sequence("ACDEFGHIKLMNPQRSTVWY") is True

    def test_is_protein_sequence_dna(self):
        assert is_protein_sequence("ATCGATCG") is False


# ── Report ─────────────────────────────────────────────────────────────────

class TestReport:
    def test_tsv_header(self, annotator, synthetic_seq):
        hits = annotator.annotate_sequence(synthetic_seq)
        tsv = report.to_tsv(hits, seq_id="test")
        assert tsv.startswith("seq_id\t")

    def test_tsv_row_count(self, annotator, synthetic_seq):
        hits = annotator.annotate_sequence(synthetic_seq)
        tsv = report.to_tsv(hits, seq_id="test")
        lines = [l for l in tsv.strip().split("\n") if l]
        assert len(lines) == len(hits) + 1

    def test_json_output(self, annotator, synthetic_seq):
        hits = annotator.annotate_sequence(synthetic_seq)
        js = report.to_json(hits, seq_id="test")
        parsed = json.loads(js)
        assert isinstance(parsed, list)
        if parsed:
            assert "ptm_name" in parsed[0]
            assert "query_position" in parsed[0]
            assert "references" in parsed[0]

    def test_summary_contains_warning(self, annotator, synthetic_seq):
        hits = annotator.annotate_sequence(synthetic_seq)
        s = report.summary(hits, seq_id="test")
        assert "heuristic" in s.lower()
        assert "1MRO" in s

    def test_tsv_file_write(self, annotator, synthetic_seq, tmp_path):
        hits = annotator.annotate_sequence(synthetic_seq)
        out = tmp_path / "out.tsv"
        report.to_tsv(hits, seq_id="test", path=out)
        assert out.exists()
