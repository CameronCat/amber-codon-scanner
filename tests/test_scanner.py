
import json
import pytest
from amber_codon_scanner import AmberCodonScanner, AmberCodon, Classification
from amber_codon_scanner.utils import reverse_complement, translate_codon, parse_fasta
from amber_codon_scanner import report


def _make_mid_orf_uag(upstream_codons=15, downstream_codons=10):
    up   = "ATG" + "GCG" * (upstream_codons - 1)
    down = "GCG" * downstream_codons + "TGA"
    return up + "TAG" + down


def _make_terminal_uag():
    return "ATG" + "GCG" * 3 + "TAG"


@pytest.fixture
def scanner():
    return AmberCodonScanner(min_orf_length=5, downstream_stop_window=20)


@pytest.fixture
def tmp_fasta(tmp_path):
    seq = _make_mid_orf_uag()
    fa = tmp_path / "test.fasta"
    fa.write_text(f">seq1\n{seq}\n")
    return fa


class TestUtils:
    def test_reverse_complement(self):
        assert reverse_complement("ATCG") == "CGAT"

    def test_translate_stop_uag(self):
        assert translate_codon("TAG") == "*"

    def test_translate_stop_taa(self):
        assert translate_codon("TAA") == "*"

    def test_translate_met(self):
        assert translate_codon("ATG") == "M"

    def test_translate_unknown(self):
        assert translate_codon("NNN") == "X"

    def test_parse_fasta(self, tmp_fasta):
        records = list(parse_fasta(tmp_fasta))
        assert len(records) == 1
        seq_id, seq = records[0]
        assert seq_id == "seq1"
        assert "TAG" in seq


class TestScanner:
    def test_finds_uag_in_sequence(self, scanner):
        seq = _make_mid_orf_uag()
        results = scanner.scan_sequence(seq, seq_id="test")
        uag_positions = [c.position for c in results if c.strand == "+"]
        assert len(uag_positions) >= 1

    def test_mid_orf_detected(self, scanner):
        seq = _make_mid_orf_uag(upstream_codons=10, downstream_codons=10)
        results = scanner.scan_sequence(seq, seq_id="test")
        fwd = [c for c in results if c.strand == "+"]
        assert any(c.mid_orf for c in fwd)

    def test_terminal_uag_not_mid_orf(self, scanner):
        seq = _make_terminal_uag()
        results = scanner.scan_sequence(seq, seq_id="test")
        fwd = [c for c in results if c.strand == "+"]
        terminal = [c for c in fwd if not c.mid_orf]
        assert len(terminal) >= 1

    def test_classification_stop_likely_for_terminal(self, scanner):
        seq = _make_terminal_uag()
        results = scanner.scan_sequence(seq, seq_id="test")
        fwd = [c for c in results if c.strand == "+"]
        stop_likely = [c for c in fwd if c.classification == Classification.STOP_LIKELY]
        assert len(stop_likely) >= 1

    def test_scan_fasta(self, scanner, tmp_fasta):
        results = scanner.scan_fasta(str(tmp_fasta))
        assert "seq1" in results
        assert len(results["seq1"]) >= 1

    def test_both_strands_scanned(self, scanner):
        seq = "ATGCTAGATG" * 6 + "TGA"
        results = scanner.scan_sequence(seq)
        strands = {c.strand for c in results}
        assert "+" in strands
        assert "-" in strands

    def test_frame_assignment(self, scanner):
        seq = "A" * 6 + "TAG" + "A" * 20 + "TGA"
        results = scanner.scan_sequence(seq)
        fwd = [c for c in results if c.strand == "+"]
        at_6 = [c for c in fwd if c.position == 6]
        assert at_6
        assert at_6[0].frame == 0

    def test_pylis_check_disabled(self):
        scanner_no_pylis = AmberCodonScanner(check_pylis=False)
        seq = _make_mid_orf_uag()
        results = scanner_no_pylis.scan_sequence(seq)
        assert all(not c.pylis_detected for c in results)

    def test_ambiguous_when_mid_orf_no_pylis(self, scanner):
        seq = "ATG" + "GCG" * 10 + "TAG" + "ATA" * 15 + "TAA"
        results = scanner.scan_sequence(seq)
        fwd = [c for c in results if c.strand == "+"]
        mid = [c for c in fwd if c.mid_orf]
        if mid:
            assert mid[0].classification in (
                Classification.AMBIGUOUS, Classification.PYL_CANDIDATE
            )


class TestClassification:
    def test_pyl_candidate_requires_both(self, scanner):
        cls, _ = scanner._classify(mid_orf=True, pylis=True, ds_gc=0.6)
        assert cls == Classification.PYL_CANDIDATE

    def test_ambiguous_pylis_no_midorf(self, scanner):
        cls, _ = scanner._classify(mid_orf=False, pylis=True, ds_gc=0.6)
        assert cls == Classification.AMBIGUOUS

    def test_ambiguous_midorf_no_pylis(self, scanner):
        cls, _ = scanner._classify(mid_orf=True, pylis=False, ds_gc=0.3)
        assert cls == Classification.AMBIGUOUS

    def test_stop_likely_neither(self, scanner):
        cls, _ = scanner._classify(mid_orf=False, pylis=False, ds_gc=0.2)
        assert cls == Classification.STOP_LIKELY


class TestReport:
    def test_tsv_header(self, scanner):
        seq = _make_mid_orf_uag()
        codons = scanner.scan_sequence(seq)
        tsv = report.to_tsv(codons)
        assert tsv.startswith("seq_id\t")

    def test_tsv_row_count(self, scanner):
        seq = _make_mid_orf_uag()
        codons = scanner.scan_sequence(seq)
        tsv = report.to_tsv(codons)
        lines = [l for l in tsv.strip().split("\n") if l]
        assert len(lines) == len(codons) + 1

    def test_json_output(self, scanner):
        seq = _make_mid_orf_uag()
        codons = scanner.scan_sequence(seq)
        js = report.to_json(codons)
        parsed = json.loads(js)
        assert isinstance(parsed, list)
        if parsed:
            assert "classification" in parsed[0]
            assert "position" in parsed[0]

    def test_summary_contains_counts(self, scanner):
        seq = _make_mid_orf_uag()
        codons = scanner.scan_sequence(seq)
        s = report.summary(codons, show_all=True)
        assert "Total UAG codons found" in s
        assert "Pyrrolysine candidates" in s

    def test_summary_warning_present(self, scanner):
        seq = _make_mid_orf_uag()
        codons = scanner.scan_sequence(seq)
        s = report.summary(codons)
        assert "heuristic" in s.lower()

    def test_tsv_file_write(self, scanner, tmp_path):
        seq = _make_mid_orf_uag()
        codons = scanner.scan_sequence(seq)
        out = tmp_path / "out.tsv"
        report.to_tsv(codons, path=out)
        assert out.exists()
        assert out.stat().st_size > 0
