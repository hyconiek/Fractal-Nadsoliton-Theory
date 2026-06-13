from __future__ import annotations

import json
import subprocess
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent
SCRIPT = ROOT / "p2666_s1616_pair12_selector_witness_to_boundary_phase_sector_descent_audit.py"
OUT = ROOT / "generated" / "p2666_s1616_pair12_selector_witness_to_boundary_phase_sector_descent_audit.json"
MD = ROOT / "generated" / "p2666_s1616_pair12_selector_witness_to_boundary_phase_sector_descent_audit.md"


class P2666Pair12SelectorWitnessToBoundaryPhaseSectorDescentAuditTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        subprocess.run([sys.executable, str(SCRIPT)], cwd=ROOT, check=True)
        cls.payload = json.loads(OUT.read_text(encoding="utf-8"))

    def test_content_first_grep_covers_pair12_selector_and_boundary_phase(self) -> None:
        audit = self.payload["semantic_rg_antiduplication_audit"]
        self.assertIn("content-first", audit["mode"])
        for key in (
            "pair12_witness_split_content",
            "typed_descent_promotion_content",
            "selector_source_content",
            "boundary_phase_sector_content",
            "nonclosure_guard_content",
        ):
            self.assertIn(key, audit["patterns"])
            self.assertGreater(audit["patterns"][key]["count"], 0)

    def test_prior_summary_audit_identifies_near_miss_and_blockers(self) -> None:
        summary = self.payload["prior_summary_audit"]
        self.assertTrue(summary["p731_w_break_separates_pair12"])
        self.assertTrue(summary["p731_pair2_positive"])
        self.assertFalse(summary["p731_promotion_bridge_exported"])
        self.assertTrue(summary["p741_actual_source_witness_exported"])
        self.assertTrue(summary["p741_witness_remains_prelm_not_pair12_typed"])
        self.assertFalse(summary["p741_promotion_bridge_exported"])
        self.assertTrue(summary["p764_typed_descent_target_exported"])
        self.assertTrue(summary["p764_typed_descent_target_future_route_only"])
        self.assertTrue(summary["p2665_no_accepted_selector_bridge"])

    def test_descent_mapping_witness_blocks_convention_false_pass(self) -> None:
        witness = self.payload["descent_mapping_witness"]
        self.assertEqual(witness["orientation_convention_degeneracy_count"], 2)
        self.assertTrue(witness["raw_holonomy_sign_preserved_by_at_least_one_mapping"])
        self.assertFalse(witness["strict_current_descent_exported_now"])
        self.assertEqual(witness["passing_descent_mappings"], [])
        rows = witness["rows"]
        self.assertEqual(len(rows), 2)
        self.assertTrue(any(row["selects_nonexact_boundary_sector_one"] for row in rows))
        self.assertTrue(all(not row["passes_descent_acceptance_now"] for row in rows))

    def test_no_closure_and_docs_updated(self) -> None:
        decision = self.payload["closure_decision"]
        self.assertTrue(decision["pair12_split_exported"])
        self.assertTrue(decision["pair2_positive_branch"])
        self.assertTrue(decision["mathematically_possible_sector_one_mapping_exists"])
        self.assertEqual(decision["passing_descent_mapping_count"], 0)
        self.assertFalse(decision["strict_current_descent_exported_now"])
        self.assertFalse(decision["boundary_phase_bit_target_exported_now"])
        self.assertFalse(decision["unconditional_uv_unit_selected_now"])
        self.assertFalse(decision["beta_source_exported_now"])
        self.assertFalse(decision["qw2191_discharged_now"])
        self.assertFalse(decision["role_bearing_ltotal_now"])
        self.assertFalse(decision["toe_closure_now"])
        self.assertTrue(all(value is False for value in self.payload["negative_export_flags"].values()))
        self.assertIn("P2666/S1616", MD.read_text(encoding="utf-8"))
        self.assertIn("P2666/S1616", (ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md").read_text(encoding="utf-8"))
        self.assertIn("P2666/S1616", (ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md").read_text(encoding="utf-8"))


if __name__ == "__main__":
    unittest.main()
