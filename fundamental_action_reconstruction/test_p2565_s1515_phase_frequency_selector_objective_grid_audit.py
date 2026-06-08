from __future__ import annotations

import json
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT))

from p2565_s1515_phase_frequency_selector_objective_grid_audit import MD, OUT, append_doc_sections, build_payload, write_markdown


class P2565PhaseFrequencySelectorObjectiveGridAuditTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.payload = build_payload()
        OUT.write_text(json.dumps(cls.payload, indent=2, sort_keys=True), encoding="utf-8")
        write_markdown(cls.payload)
        append_doc_sections()
        cls.theorem = cls.payload["phase_frequency_selector_objective_grid_audit"]["theorem_export"]
        cls.scan = cls.theorem["selector_objective_grid_audit"]

    def test_identity_and_inheritance(self) -> None:
        self.assertEqual(self.payload["packet_id"], "P2565")
        self.assertEqual(self.payload["stage_id"], "S1515")
        self.assertIn("SELECTOR_OBJECTIVE_GRID_AUDIT", self.payload["status"])
        self.assertEqual(self.theorem["frontier_atom_under_attack"], "strict_phase_frequency_source")
        self.assertTrue(self.theorem["p2564_open_sign_cell_inherited"])
        self.assertTrue(self.theorem["p2563_rational_winding_obstruction_inherited"])
        self.assertTrue(self.theorem["p2561_phase_frequency_residual_atom_inherited"])

    def test_objective_grid_results(self) -> None:
        self.assertEqual(self.scan["grid_steps_per_axis"], 21)
        self.assertEqual(self.scan["candidate_grid_point_count"], 441)
        self.assertEqual(self.scan["accepted_sign_preserving_point_count"], 441)
        self.assertEqual(self.scan["objective_count"], 5)
        self.assertEqual(self.scan["objectives_with_non_strict_winner_count"], 5)
        self.assertTrue(self.theorem["all_audited_objectives_have_non_strict_winners"])
        self.assertTrue(self.theorem["source_free_objective_choice_remains_extra_obligation"])
        for winner in self.scan["objective_winners"].values():
            self.assertFalse(winner["is_strict_tuple"])
            self.assertGreater(winner["improvement_over_strict"], 0.0)

    def test_toe_and_closure_diagnosis(self) -> None:
        self.assertIn("positive_potential", self.theorem["toe_potential_assessment"])
        self.assertIn("current_limit", self.theorem["toe_potential_assessment"])
        self.assertIn("why_closure_is_hard", self.theorem["closure_problem_diagnosis"])
        self.assertGreaterEqual(len(self.theorem["closure_problem_diagnosis"]["why_closure_is_hard"]), 5)
        self.assertIn("derive one selector functional", self.theorem["recommended_next_honest_step"])

    def test_negative_controls_and_docs(self) -> None:
        self.assertFalse(self.theorem["selector_objective_source_exported"])
        self.assertFalse(self.theorem["strict_phase_frequency_source_exported"])
        self.assertFalse(self.theorem["strict_dynamical_source_for_A_P_D_exported"])
        self.assertFalse(self.theorem["strict_damping_beta_eta_source_exported"])
        self.assertFalse(self.theorem["bridge_theorem_exported"])
        self.assertFalse(self.theorem["role_transfer_theorem_exported"])
        self.assertFalse(self.theorem["qw2191_discharged_by_this_certificate"])
        self.assertFalse(self.theorem["toe_closure_claimed"])
        self.assertTrue(all(self.payload["gatekeeper_checks"].values()))
        self.assertEqual(self.payload["rg_non_duplication_audit"]["tool"], "rg")
        self.assertIn("P2565/S1515", MD.read_text(encoding="utf-8"))
        self.assertIn("P2565/S1515", (ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md").read_text(encoding="utf-8"))
        self.assertIn("P2565/S1515", (ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md").read_text(encoding="utf-8"))


if __name__ == "__main__":
    unittest.main()
