from __future__ import annotations

import json
import subprocess
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent
SCRIPT = ROOT / "p2669_s1619_boundary_cycle_cross_invariant_boolean_ansatz_audit.py"
OUT = ROOT / "generated" / "p2669_s1619_boundary_cycle_cross_invariant_boolean_ansatz_audit.json"
MD = ROOT / "generated" / "p2669_s1619_boundary_cycle_cross_invariant_boolean_ansatz_audit.md"


class P2669BoundaryCycleCrossInvariantBooleanAnsatzAuditTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        subprocess.run([sys.executable, str(SCRIPT)], cwd=ROOT, check=True)
        cls.payload = json.loads(OUT.read_text(encoding="utf-8"))

    def test_content_first_grep_covers_boolean_cross_invariant(self) -> None:
        audit = self.payload["semantic_rg_antiduplication_audit"]
        self.assertIn("content-first", audit["mode"])
        for key in (
            "cross_invariant_boolean_content",
            "boundary_cycle_sector_content",
            "pair2_selector_content",
            "source_origin_guard_content",
            "nonclosure_guard_content",
        ):
            self.assertIn(key, audit["patterns"])
            self.assertGreater(audit["patterns"][key]["count"], 0)

    def test_upstream_p2668_consistency(self) -> None:
        upstream = self.payload["upstream_consistency"]
        self.assertTrue(upstream["p2668_existing_orientation_export_present"])
        self.assertTrue(upstream["p2668_no_source_forbids_sector_swap"])
        self.assertTrue(upstream["p2668_no_source_ties_pair2_to_sector1"])
        self.assertTrue(upstream["p2668_no_boundary_phase_bit_target"])

    def test_boolean_ansatz_finds_mathematical_candidates_but_no_source(self) -> None:
        witness = self.payload["boolean_ansatz_witness"]
        self.assertEqual(witness["total_boolean_ansatz_count"], 16)
        self.assertTrue(witness["mathematical_candidates_exist"])
        self.assertFalse(witness["mixed_candidates_exist"])
        self.assertGreater(witness["mathematical_sector_swap_odd_tie_count"], 0)
        self.assertEqual(witness["mixed_pair_sector_candidate_count"], 0)
        self.assertEqual(witness["passing_cross_invariant_count"], 0)
        self.assertFalse(witness["convention_free_source_exported_for_any_candidate"])
        self.assertTrue(any(row["distinguishes_pair_branch"] for row in witness["rows"] if row["sector_swap_odd"] and row["ties_pair2_positive_to_sector1"]))

    def test_no_closure_and_docs_updated(self) -> None:
        decision = self.payload["closure_decision"]
        self.assertTrue(decision["mathematical_candidates_exist"])
        self.assertFalse(decision["mixed_candidates_exist"])
        self.assertEqual(decision["passing_cross_invariant_count"], 0)
        self.assertFalse(decision["convention_free_source_exported_for_any_candidate"])
        self.assertFalse(decision["boundary_phase_bit_target_exported_now"])
        self.assertFalse(decision["unconditional_uv_unit_selected_now"])
        self.assertFalse(decision["beta_source_exported_now"])
        self.assertFalse(decision["qw2191_discharged_now"])
        self.assertFalse(decision["role_bearing_ltotal_now"])
        self.assertFalse(decision["toe_closure_now"])
        self.assertTrue(all(value is False for value in self.payload["negative_export_flags"].values()))
        self.assertIn("P2669/S1619", MD.read_text(encoding="utf-8"))
        self.assertIn("P2669/S1619", (ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md").read_text(encoding="utf-8"))
        self.assertIn("P2669/S1619", (ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md").read_text(encoding="utf-8"))


if __name__ == "__main__":
    unittest.main()
