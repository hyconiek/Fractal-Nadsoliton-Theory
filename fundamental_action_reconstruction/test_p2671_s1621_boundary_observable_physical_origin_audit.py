from __future__ import annotations

import json
import subprocess
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent
SCRIPT = ROOT / "p2671_s1621_boundary_observable_physical_origin_audit.py"
OUT = ROOT / "generated" / "p2671_s1621_boundary_observable_physical_origin_audit.json"
MD = ROOT / "generated" / "p2671_s1621_boundary_observable_physical_origin_audit.md"


class P2671BoundaryObservablePhysicalOriginAuditTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        subprocess.run([sys.executable, str(SCRIPT)], cwd=ROOT, check=True)
        cls.payload = json.loads(OUT.read_text(encoding="utf-8"))

    def test_content_first_grep_covers_observable_origin_and_selector(self) -> None:
        audit = self.payload["semantic_rg_antiduplication_audit"]
        self.assertIn("content-first", audit["mode"])
        for key in (
            "observable_origin_content",
            "pair_variable_content",
            "boundary_sector_content",
            "auxiliary_lift_content",
            "selector_content",
            "nonclosure_guard_content",
        ):
            self.assertIn(key, audit["patterns"])
            self.assertGreater(audit["patterns"][key]["count"], 0)

    def test_upstream_p2670_consistency(self) -> None:
        upstream = self.payload["upstream_consistency"]
        self.assertTrue(upstream["p2670_higher_order_candidates_exist"])
        self.assertTrue(upstream["p2670_no_passing_source"])
        self.assertTrue(upstream["p2670_no_boundary_phase_bit_target"])
        self.assertTrue(upstream["p2670_no_ltotal_reopening"])

    def test_candidate_matrix_finds_near_misses_but_no_single_origin(self) -> None:
        matrix = self.payload["observable_candidate_matrix"]
        self.assertGreaterEqual(len(matrix), 6)
        self.assertTrue(any(row["defines_pair_variable"] for row in matrix))
        self.assertTrue(any(row["defines_boundary_sector_variable"] for row in matrix))
        self.assertTrue(any(row["defines_auxiliary_lift"] for row in matrix))
        self.assertTrue(all(row["content_hits"] > 0 for row in matrix))
        self.assertFalse(any(row["passes_single_observable_origin_now"] for row in matrix))

    def test_subset_witness_blocks_combined_origin_false_pass(self) -> None:
        witness = self.payload["subset_origin_witness"]
        self.assertEqual(witness["total_subsets_checked"], 63)
        self.assertGreater(witness["near_miss_pair_and_sector_subset_count"], 0)
        self.assertEqual(witness["passing_subset_count"], 0)
        self.assertEqual(witness["passing_subsets"], [])

    def test_no_closure_and_docs_updated(self) -> None:
        decision = self.payload["closure_decision"]
        self.assertEqual(decision["passing_subset_count"], 0)
        self.assertGreater(decision["near_miss_pair_and_sector_subset_count"], 0)
        self.assertFalse(decision["boundary_phase_bit_target_exported_now"])
        self.assertFalse(decision["beta_source_exported_now"])
        self.assertFalse(decision["qw2191_discharged_now"])
        self.assertFalse(decision["role_bearing_ltotal_now"])
        self.assertFalse(decision["toe_closure_now"])
        self.assertTrue(all(value is False for value in self.payload["negative_export_flags"].values()))
        self.assertIn("P2671/S1621", MD.read_text(encoding="utf-8"))
        self.assertIn("P2671/S1621", (ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md").read_text(encoding="utf-8"))
        self.assertIn("P2671/S1621", (ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md").read_text(encoding="utf-8"))


if __name__ == "__main__":
    unittest.main()
