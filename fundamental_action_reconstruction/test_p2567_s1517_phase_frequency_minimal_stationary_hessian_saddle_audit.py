from __future__ import annotations

import json
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT))

from p2567_s1517_phase_frequency_minimal_stationary_hessian_saddle_audit import MD, OUT, append_doc_sections, build_payload, write_markdown


class P2567PhaseFrequencyMinimalStationaryHessianSaddleAuditTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.payload = build_payload()
        OUT.write_text(json.dumps(cls.payload, indent=2, sort_keys=True), encoding="utf-8")
        write_markdown(cls.payload)
        append_doc_sections()
        cls.theorem = cls.payload["phase_frequency_minimal_stationary_hessian_saddle_audit"]["theorem_export"]
        cls.audit = cls.theorem["minimal_three_node_stationary_hessian_audit"]

    def test_identity_and_inheritance(self) -> None:
        self.assertEqual(self.payload["packet_id"], "P2567")
        self.assertEqual(self.payload["stage_id"], "S1517")
        self.assertIn("MINIMAL_STATIONARY_HESSIAN_SADDLE", self.payload["status"])
        self.assertEqual(self.theorem["frontier_atom_under_attack"], "strict_phase_frequency_source")
        self.assertTrue(self.theorem["p2566_stationarity_underidentification_inherited"])
        self.assertTrue(self.theorem["p2565_objective_choice_obligation_inherited"])
        self.assertTrue(self.theorem["p2561_phase_frequency_residual_atom_inherited"])

    def test_all_minimal_stationary_witnesses_are_saddles(self) -> None:
        self.assertEqual(self.audit["support_size"], 3)
        self.assertEqual(self.audit["audited_support_count"], 220)
        self.assertEqual(self.audit["expected_support_count_choose_12_3"], 220)
        self.assertEqual(self.audit["indefinite_saddle_count"], 220)
        self.assertEqual(self.audit["non_saddle_count"], 0)
        self.assertTrue(self.audit["all_minimal_three_node_stationary_witnesses_are_saddles"])
        self.assertLess(self.audit["max_stationarity_residual_abs"], 1e-12)
        self.assertTrue(self.theorem["minimal_signed_stationary_witnesses_fail_second_order_selector_test"])
        for row in self.audit["sample_saddle_rows"]:
            self.assertLess(row["determinant"], 0.0)
            self.assertEqual(row["classification"], "indefinite_saddle")

    def test_negative_controls_and_docs(self) -> None:
        self.assertFalse(self.theorem["minimal_stationary_hessian_selector_source_exported"])
        self.assertFalse(self.theorem["strict_phase_frequency_source_exported"])
        self.assertFalse(self.theorem["bridge_theorem_exported"])
        self.assertFalse(self.theorem["role_transfer_theorem_exported"])
        self.assertFalse(self.theorem["qw2191_discharged_by_this_certificate"])
        self.assertFalse(self.theorem["toe_closure_claimed"])
        self.assertTrue(all(self.payload["gatekeeper_checks"].values()))
        self.assertEqual(self.payload["rg_non_duplication_audit"]["tool"], "rg")
        self.assertIn("P2567/S1517", MD.read_text(encoding="utf-8"))
        self.assertIn("P2567/S1517", (ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md").read_text(encoding="utf-8"))
        self.assertIn("P2567/S1517", (ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md").read_text(encoding="utf-8"))


if __name__ == "__main__":
    unittest.main()
