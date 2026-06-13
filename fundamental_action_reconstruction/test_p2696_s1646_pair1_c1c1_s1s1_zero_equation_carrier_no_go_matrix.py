from __future__ import annotations

import json
import subprocess
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent
SCRIPT = ROOT / "p2696_s1646_pair1_c1c1_s1s1_zero_equation_carrier_no_go_matrix.py"
OUT = ROOT / "generated" / "p2696_s1646_pair1_c1c1_s1s1_zero_equation_carrier_no_go_matrix.json"
MD = ROOT / "generated" / "p2696_s1646_pair1_c1c1_s1s1_zero_equation_carrier_no_go_matrix.md"


class P2696Pair1C1C1S1S1ZeroEquationCarrierNoGoMatrixTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        subprocess.run([sys.executable, str(SCRIPT)], cwd=ROOT, check=True)
        cls.payload = json.loads(OUT.read_text(encoding="utf-8"))

    def test_content_grep_covers_pair1_no_go_inputs(self) -> None:
        audit = self.payload["content_grep"]
        self.assertIn("pair1 c1c1/s1s1", audit["mode"])
        for key in ("p2695_selected_p2696", "value_instance_obstruction", "finite_sign_scan_obstruction", "p631_negative_freeze", "forbidden_replays"):
            self.assertIn(key, audit["patterns"])
            self.assertGreater(audit["patterns"][key]["count"], 0)

    def test_state_reads_current_obstructions(self) -> None:
        state = self.payload["state_reads"]
        self.assertTrue(state["p2695_recommended_p2696"])
        self.assertTrue(state["p2695_remaining_mentions_c1c1_s1s1"])
        self.assertEqual(set(state["p477_violated_equations"]), {"c1c1", "s1s1"})
        self.assertFalse(state["p477_all_zero_equations_satisfied"])
        self.assertTrue(state["n520_states_value_instance_obstruction"])
        self.assertFalse(state["p479_any_reference_has_solution"])
        self.assertTrue(state["n522_finite_family_obstruction_discharged"])
        self.assertTrue(state["p631_negative_freeze_selected"])

    def test_zero_equation_matrix_is_two_row_bounded_no_go(self) -> None:
        rows = self.payload["zero_equation_matrix"]
        self.assertEqual(len(rows), 2)
        self.assertEqual({row["target_id"] for row in rows}, {"declared_pair1_c1c1_zero_equation_carrier", "declared_pair1_s1s1_zero_equation_carrier"})
        self.assertTrue(all(row["p477_marks_violated"] for row in rows))
        self.assertTrue(all(not row["p477_is_zero"] for row in rows))
        self.assertTrue(all(row["n520_value_instance_obstruction_supports_no_go"] for row in rows))
        self.assertTrue(all(row["p479_n522_finite_sign_scan_no_solution_supports_no_go"] for row in rows))
        self.assertTrue(all(not row["explicit_zero_witness_exported_now"] for row in rows))

    def test_lane_boundary_and_next_step(self) -> None:
        boundary = self.payload["lane_boundary"]
        self.assertTrue(boundary["all_targeted_pair1_c1c1_s1s1_rows_bounded_no_go"])
        self.assertFalse(boundary["strict_zero_witness_exported_for_c1c1_or_s1s1"])
        self.assertTrue(boundary["finite_sign_scan_blocks_reference_family_only"])
        self.assertTrue(boundary["direct_formal_residual_cancellation_lane_frozen_by_p631"])
        self.assertTrue(boundary["qw2191_still_separate_open_boundary"])
        self.assertFalse(boundary["full_route_closure_exported"])
        self.assertFalse(boundary["strict_core_promotion_exported"])
        decision = self.payload["decision"]
        self.assertIn("P2697", decision["next_honest_step"])
        self.assertFalse(decision["direct_route_reopen_admissible_without_new_provider"])
        self.assertTrue(all(value is False for value in self.payload["negative_export_flags"].values()))
        self.assertIn("P2696/S1646", MD.read_text(encoding="utf-8"))
        self.assertIn("P2696/S1646", (ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md").read_text(encoding="utf-8"))
        self.assertIn("P2696/S1646", (ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md").read_text(encoding="utf-8"))
        self.assertIn("P2696/S1646", (ROOT.parent / "AGENTS.md").read_text(encoding="utf-8"))


if __name__ == "__main__":
    unittest.main()
