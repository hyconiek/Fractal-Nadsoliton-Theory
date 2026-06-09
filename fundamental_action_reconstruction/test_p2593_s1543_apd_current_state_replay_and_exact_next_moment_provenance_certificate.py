from __future__ import annotations

import json
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT))

from p2593_s1543_apd_current_state_replay_and_exact_next_moment_provenance_certificate import MD, OUT, append_doc_sections, build_payload, write_markdown


class P2593APDCurrentStateReplayExactNextMomentProvenanceCertificateTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.payload = build_payload()
        OUT.write_text(json.dumps(cls.payload, indent=2, sort_keys=True), encoding="utf-8")
        write_markdown(cls.payload)
        append_doc_sections()
        cls.theorem = cls.payload["apd_current_state_replay_and_exact_next_moment_provenance_certificate"]["theorem_export"]
        cls.exact = cls.theorem["exact_newton_girard_replay"]

    def test_identity_current_state_and_inheritance(self) -> None:
        self.assertEqual(self.payload["packet_id"], "P2593")
        self.assertEqual(self.payload["stage_id"], "S1543")
        self.assertIn("CURRENT_STATE_REPLAY", self.payload["status"])
        self.assertEqual(self.theorem["frontier_atom_under_attack"], "strict_dynamical_source_for_A_P_D")
        self.assertTrue(self.theorem["p2592_next_even_moment_certificate_inherited"])
        self.assertTrue(self.theorem["p2592_whole_source_artifact_mismatch_detected"] or self.theorem["p2592_whole_source_artifact_matches_current_repo"])
        self.assertTrue(self.theorem["p2592_theorem_relevant_inputs_replay_on_current_repo"])
        self.assertTrue(self.theorem["p2592_commit_is_on_current_head_history"])

    def test_exact_replay_refines_p2592(self) -> None:
        self.assertEqual(self.exact["newton_girard_internal_formula_exact"], "74658 - 4*e4")
        self.assertEqual(self.exact["central_eighth_formula_exact_rational"], "42715646049/32768 - 8*e4")
        self.assertEqual(self.exact["recorded_internal_formula"], "74658.0 - 4*e4")
        self.assertTrue(self.exact["internal_formula_replays_recorded"])
        self.assertTrue(self.exact["central_float_formula_replays_recorded"])
        self.assertEqual(self.exact["max_abs_internal_grid_delta_vs_recorded"], 0.0)
        self.assertLessEqual(self.exact["max_abs_central_grid_delta_vs_recorded"], 1.0e-9)
        self.assertTrue(self.theorem["exact_replay_preserves_p2592_numeric_certificate"])
        stability = self.exact["recomputed_gatekeeper_stability"]
        self.assertTrue(stability["lower_shells_constant_but_next_shell_varies"])
        self.assertEqual(stability["internal_eighth_shell_slope_by_product_parameter"], -4.0)
        self.assertEqual(stability["central_eighth_moment_slope_by_product_parameter"], -8.0)

    def test_negative_controls_and_docs(self) -> None:
        self.assertFalse(self.theorem["apd_current_state_replay_source_exported"])
        self.assertFalse(self.theorem["apd_provenance_replay_selector_source_exported"])
        self.assertFalse(self.theorem["apd_exact_rational_next_moment_source_exported"])
        self.assertFalse(self.theorem["strict_dynamical_source_for_A_P_D_exported"])
        self.assertFalse(self.theorem["bridge_theorem_exported"])
        self.assertFalse(self.theorem["role_transfer_theorem_exported"])
        self.assertFalse(self.theorem["qw2191_discharged_by_this_certificate"])
        self.assertFalse(self.theorem["toe_closure_claimed"])
        self.assertTrue(all(self.payload["gatekeeper_checks"].values()))
        self.assertEqual(self.payload["rg_non_duplication_audit"]["tool"], "rg")
        self.assertIn("P2593/S1543", MD.read_text(encoding="utf-8"))
        self.assertIn("P2593/S1543", (ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md").read_text(encoding="utf-8"))
        self.assertIn("P2593/S1543", (ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md").read_text(encoding="utf-8"))


if __name__ == "__main__":
    unittest.main()
