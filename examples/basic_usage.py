
from amber_codon_scanner import AmberCodonScanner, Classification, report

print("=" * 60)
print("Example 1: Mid-ORF UAG (pyrrolysine candidate context)")
print("=" * 60)

mid_orf_seq = "ATG" + "GCG" * 15 + "TAG" + "GCG" * 12 + "TGA"

scanner = AmberCodonScanner(min_orf_length=5, downstream_stop_window=20)
codons = scanner.scan_sequence(mid_orf_seq, seq_id="synthetic_mid_orf")
print(report.summary(codons, show_all=True))

print("\n" + "=" * 60)
print("Example 2: Terminal UAG (stop codon context)")
print("=" * 60)

terminal_seq = "ATG" + "GCG" * 5 + "TAG"
codons_term = scanner.scan_sequence(terminal_seq, seq_id="synthetic_terminal")
print(report.summary(codons_term, show_all=True))

print("\n" + "=" * 60)
print("Example 3: Scan from FASTA file")
print("=" * 60)
results = scanner.scan_fasta("examples/example_sequences.fasta")
for seq_id, codons in results.items():
    print(f"\n--- {seq_id} ---")
    print(report.summary(codons, show_all=True))

print("\n" + "=" * 60)
print("Example 4: Export results")
print("=" * 60)

import json
all_codons = []
for c_list in results.values():
    all_codons.extend(c_list)

tsv = report.to_tsv(all_codons)
print("TSV (first 3 lines):")
for line in tsv.split("\n")[:3]:
    print(" ", line)

js = report.to_json(all_codons)
parsed = json.loads(js)
print(f"\nJSON: {len(parsed)} UAG codons serialised")

print("\n" + "=" * 60)
print("Example 5: Filter to pyrrolysine candidates only")
print("=" * 60)

candidates = [c for c in all_codons if c.classification == Classification.PYL_CANDIDATE]
print(f"Found {len(candidates)} pyrrolysine candidate(s)")
for c in candidates:
    print(f"  {c.seq_id}  pos={c.position}  strand={c.strand}  "
          f"PYLIS={c.pylis_detected}  mid_orf={c.mid_orf}")

print("\nDone.")
print("\n⚠ Remember: classification is heuristic only.")
print("  Confirm with pylTSBCD gene cluster presence and/or experimental data.")
