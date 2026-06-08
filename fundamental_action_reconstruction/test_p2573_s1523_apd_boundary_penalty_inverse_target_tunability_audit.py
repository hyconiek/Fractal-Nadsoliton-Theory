from __future__ import annotations

import json
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT))

from p2573_s1523_apd_boundary_penalty_inverse_target_tunability_audit import MD, OUT, append_doc_sections, build_payload, write_markdown


class P2573APDBoundaryPenaltyInverseTargetTunabilityTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.payload = build_payload()
        OUT.write_text(json.dumps(cls.payload, indent=2, sort_keys=True), encoding="utf-8")
        write_markdown(cls.payload)
        append_doc_sections()
        cls.theorem = cls.payload["apd_boundary_penalty_inverse_target_tunability_audit"]["theorem_export"]
        cls.audit = cls.theorem["boundary_penalty_inverse_target_tunability_audit"]

    def test_identity_and_inheritance(self) -> None:
        self.assertEqual(self.payload["packet_id"], "P2573")
        self.assertEqual(self.payload["stage_id"], "S1523")
        self.assertIn("BOUNDARY_PENALTY_INVERSE_TARGET_TUNABILITY", self.payload["status"])
        self.assertEqual(self.theorem["frontier_atom_under_attack"], "strict_dynamical_source_for_A_P_D")
        self.assertTrue(self.theorem["p2416_apd_value_bridge_inherited"])
        self.assertTrue(self.theorem["p2571_measure_boundary_obligation_inherited"])
        self.assertTrue(self.theorem["p2572_boundary_penalty_continuum_inherited"])

    def test_inverse_target_tunability(self) -> None:
        self.assertEqual(len(self.theorem["apd_node_rows"]), 12)
        self.assertEqual(self.audit["inverse_rows_count"], 10)
        self.assertIn("t(lambda_target)", self.audit["inverse_stationarity_formula"])
        self.assertTrue(self.audit["all_targets_recovered_with_small_error"])
        self.assertTrue(self.audit["all_solved_penalties_nonnegative"])
        self.assertTrue(self.audit["all_inverse_rows_preserve_apd_nodes"])
        self.assertTrue(self.audit["inverse_boundary_penalty_tunes_selector_post_hoc"])
        self.assertTrue(self.theorem["finite_apd_values_do_not_select_target_or_penalty_law"])
        self.assertTrue(self.theorem["inverse_penalty_fit_is_not_source_theorem"])

    def test_negative_controls_and_docs(self) -> None:
        self.assertFalse(self.theorem["apd_boundary_penalty_law_source_exported"])
        self.assertFalse(self.theorem["apd_inverse_boundary_selector_source_exported"])
        self.assertFalse(self.theorem["apd_target_lambda_source_exported"])
        self.assertFalse(self.theorem["strict_dynamical_source_for_A_P_D_exported"])
        self.assertFalse(self.theorem["bridge_theorem_exported"])
        self.assertFalse(self.theorem["role_transfer_theorem_exported"])
        self.assertFalse(self.theorem["qw2191_discharged_by_this_certificate"])
        self.assertFalse(self.theorem["toe_closure_claimed"])
        self.assertTrue(all(self.payload["gatekeeper_checks"].values()))
        self.assertEqual(self.payload["rg_non_duplication_audit"]["tool"], "rg")
        self.assertIn("P2573/S1523", MD.read_text(encoding="utf-8"))
        self.assertIn("P2573/S1523", (ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md").read_text(encoding="utf-8"))
        self.assertIn("P2573/S1523", (ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md").read_text(encoding="utf-8"))


if __name__ == "__main__":
    unittest.main()
