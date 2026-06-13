from __future__ import annotations

import json
import subprocess
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent
SCRIPT = ROOT / "p2695_s1645_residual_direct_g_family_c1s1_zero_witness_no_go_matrix.py"
OUT = ROOT / "generated" / "p2695_s1645_residual_direct_g_family_c1s1_zero_witness_no_go_matrix.json"
MD = ROOT / "generated" / "p2695_s1645_residual_direct_g_family_c1s1_zero_witness_no_go_matrix.md"


class P2695ResidualDirectGFamilyC1S1ZeroWitnessMatrixTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        subprocess.run([sys.executable, str(SCRIPT)], cwd=ROOT, check=True)
        cls.payload = json.loads(OUT.read_text(encoding="utf-8"))

    def test_content_grep_covers_selected_matrix(self) -> None:
        audit = self.payload["content_grep"]
        self.assertIn("residual direct g-family", audit["mode"])
        for key in ("p2694_selected_p2695", "g4_g6_prior_witnesses", "gy_prior_witness", "remaining_route_blockers", "forbidden_reopen"):
            self.assertIn(key, audit["patterns"])
            self.assertGreater(audit["patterns"][key]["count"], 0)

    def test_state_reads_prior_witnesses(self) -> None:
        state = self.payload["state_reads"]
        self.assertTrue(state["p2694_selected_p2695"])
        self.assertEqual(len(state["p2694_targets"]), 3)
        self.assertTrue(state["r82_g4_zero_witness"])
        self.assertTrue(state["r82_g6_zero_witness"])
        self.assertTrue(state["r83_gy_zero_witness_under_elimination"])
        self.assertTrue(state["p629_integrates_g4_g6"])
        self.assertTrue(state["p630_integrates_gy"])
        self.assertFalse(state["p630_full_closure_pass"])
        self.assertFalse(state["p630_strict_core_promotion"])

    def test_decimal_gy_witness_is_exact_zero(self) -> None:
        gy = self.payload["gy_decimal_zero_witness_check"]
        self.assertEqual(gy["shift_defect_value_decimal"], "0")
        self.assertTrue(gy["is_exact_zero_in_exported_decimal_string"])
        self.assertIn("Yukawa-elimination", gy["scope"])

    def test_residual_matrix_closes_exactly_three_g_targets(self) -> None:
        rows = self.payload["residual_g_family_matrix"]
        self.assertEqual(len(rows), 3)
        self.assertEqual({row["target_id"] for row in rows}, {"direct_g4_c1s1_shift_defect_zero_witness", "direct_g6_c1s1_shift_defect_zero_witness", "direct_gY_c1s1_shift_defect_zero_witness"})
        self.assertTrue(all(row["explicit_witness_present"] for row in rows))
        self.assertTrue(all(row["integrated_route_state"] for row in rows))
        self.assertTrue(any(row["target_id"] == "direct_gY_c1s1_shift_defect_zero_witness" and row["shift_defect_value"] == "0" for row in rows))

    def test_boundary_and_next_step_no_false_promotion(self) -> None:
        boundary = self.payload["route_boundary"]
        self.assertTrue(boundary["all_three_selected_g_family_targets_closed"])
        self.assertTrue(boundary["pair1_c1c1_zero_equation_still_unproved"])
        self.assertTrue(boundary["pair1_s1s1_zero_equation_still_unproved"])
        self.assertTrue(boundary["qw2191_selector_boundary_still_open"])
        self.assertFalse(boundary["full_route_closure_exported"])
        self.assertFalse(boundary["strict_core_promotion_exported"])
        decision = self.payload["decision"]
        self.assertIn("P2696", decision["next_honest_step"])
        self.assertTrue(decision["bounded_no_go_for_full_direct_route_now"])
        self.assertTrue(all(value is False for value in self.payload["negative_export_flags"].values()))
        self.assertIn("P2695/S1645", MD.read_text(encoding="utf-8"))
        self.assertIn("P2695/S1645", (ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md").read_text(encoding="utf-8"))
        self.assertIn("P2695/S1645", (ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md").read_text(encoding="utf-8"))
        self.assertIn("P2695/S1645", (ROOT.parent / "AGENTS.md").read_text(encoding="utf-8"))


if __name__ == "__main__":
    unittest.main()
