from __future__ import annotations

import json
import subprocess
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent
SCRIPT = ROOT / "p2646_s1596_frozen_kernel_compression_signature_preregistration_certificate.py"
OUT = ROOT / "generated" / "p2646_s1596_frozen_kernel_compression_signature_preregistration_certificate.json"
MD = ROOT / "generated" / "p2646_s1596_frozen_kernel_compression_signature_preregistration_certificate.md"


class P2646FrozenKernelCompressionSignaturePreregistrationTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        subprocess.run([sys.executable, str(SCRIPT)], cwd=ROOT, check=True)
        cls.payload = json.loads(OUT.read_text(encoding="utf-8"))

    def test_content_first_grep_and_script_correction_note_are_present(self) -> None:
        audit = self.payload["semantic_rg_antiduplication_audit"]
        self.assertIn("content-first", audit["mode"])
        self.assertIn("blind_frozen_empirical_content", audit["patterns"])
        self.assertIn("compression_signature_content", audit["patterns"])
        self.assertTrue(all(data["count"] > 0 for data in audit["patterns"].values()))
        self.assertIn("P2646/S1596", self.payload["script_correction_note"])

    def test_frozen_discriminator_gaps_are_locked_and_positive(self) -> None:
        tests = self.payload["preregistered_frozen_kernel_tests"]
        self.assertTrue(tests["all_ratio_gaps_positive"])
        self.assertTrue(tests["all_slope_gaps_positive"])
        self.assertGreater(tests["minimum_ratio_gap"], 0.0)
        self.assertGreater(tests["minimum_slope_gap"], 0.0)
        row_1_7 = next(row for row in tests["rows"] if row["near"] == 1 and row["far"] == 7)
        self.assertLess(row_1_7["strict_denominator_tail_ratio"], row_1_7["midpoint_ratio_threshold"])
        self.assertGreater(row_1_7["strict_log_tail_slope"], row_1_7["midpoint_slope_threshold"])
        self.assertLess(row_1_7["strict_over_legacy_tail_ratio"], 0.2)

    def test_preregistration_is_not_empirical_or_source_closure(self) -> None:
        decision = self.payload["closure_decision"]
        self.assertIn("PREREGISTERED_FROZEN_COMPRESSION_SIGNATURE", decision["decision"])
        self.assertFalse(decision["full_kernel_now"])
        self.assertFalse(decision["toe_closure_now"])
        self.assertTrue(decision["can_update_role_transfer_table"])
        self.assertTrue(all(value is False for value in self.payload["negative_export_flags"].values()))

    def test_upstream_and_docs_are_updated(self) -> None:
        upstream = self.payload["upstream_consistency"]
        self.assertTrue(upstream["p2644_modified_successor_not_full_transfer"])
        self.assertTrue(upstream["p2643_threshold_rejection_available"])
        self.assertTrue(upstream["p2635_empirical_interface_available"])
        self.assertTrue(upstream["p2645_role_matrix_points_to_holdout"])
        self.assertIn("P2646/S1596", MD.read_text(encoding="utf-8"))
        self.assertIn("P2646/S1596", (ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md").read_text(encoding="utf-8"))
        self.assertIn("P2646/S1596", (ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md").read_text(encoding="utf-8"))


if __name__ == "__main__":
    unittest.main()
