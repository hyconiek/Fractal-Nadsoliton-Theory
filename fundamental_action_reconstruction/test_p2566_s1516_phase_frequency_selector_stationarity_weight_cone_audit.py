from __future__ import annotations

import json
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT))

from p2566_s1516_phase_frequency_selector_stationarity_weight_cone_audit import MD, OUT, append_doc_sections, build_payload, write_markdown


class P2566PhaseFrequencySelectorStationarityWeightConeAuditTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.payload = build_payload()
        OUT.write_text(json.dumps(cls.payload, indent=2, sort_keys=True), encoding="utf-8")
        write_markdown(cls.payload)
        append_doc_sections()
        cls.theorem = cls.payload["phase_frequency_selector_stationarity_weight_cone_audit"]["theorem_export"]

    def test_identity_and_inheritance(self) -> None:
        self.assertEqual(self.payload["packet_id"], "P2566")
        self.assertEqual(self.payload["stage_id"], "S1516")
        self.assertIn("SELECTOR_STATIONARITY_WEIGHT_CONE", self.payload["status"])
        self.assertEqual(self.theorem["frontier_atom_under_attack"], "strict_phase_frequency_source")
        self.assertTrue(self.theorem["p2565_objective_choice_obligation_inherited"])
        self.assertTrue(self.theorem["p2564_open_sign_cell_inherited"])
        self.assertTrue(self.theorem["p2561_phase_frequency_residual_atom_inherited"])

    def test_rank_nullity_and_cone_obstruction(self) -> None:
        linear = self.theorem["linear_stationarity_rank_audit"]
        self.assertEqual(linear["stationarity_matrix_shape"], [2, 12])
        self.assertEqual(linear["stationarity_matrix_rank"], 2)
        self.assertEqual(linear["stationarity_nullity"], 10)
        cone = self.theorem["nonnegative_weight_cone_obstruction"]
        self.assertTrue(cone["all_sin_theta_positive_on_domain"])
        self.assertEqual(cone["positive_sign_d_range"], [0, 7])
        self.assertEqual(cone["negative_sign_d_range"], [8, 11])
        self.assertFalse(cone["nonzero_nonnegative_weight_stationarity_possible"])
        self.assertTrue(self.theorem["nonnegative_natural_weight_stationarity_rejected"])

    def test_signed_witness_and_negative_controls(self) -> None:
        witness = self.theorem["signed_weight_stationarity_witness"]
        self.assertEqual(witness["support"], [0, 7, 8])
        self.assertTrue(witness["has_negative_weight"])
        self.assertLess(witness["gradient_phi_residual_abs"], 1e-12)
        self.assertLess(witness["gradient_omega_residual_abs"], 1e-12)
        self.assertTrue(self.theorem["first_order_stationarity_with_unconstrained_signed_weights_underidentified"])
        self.assertFalse(self.theorem["stationarity_alone_selects_unique_source"])
        self.assertFalse(self.theorem["stationarity_weight_source_exported"])
        self.assertFalse(self.theorem["strict_phase_frequency_source_exported"])
        self.assertFalse(self.theorem["bridge_theorem_exported"])
        self.assertFalse(self.theorem["role_transfer_theorem_exported"])
        self.assertFalse(self.theorem["qw2191_discharged_by_this_certificate"])
        self.assertFalse(self.theorem["toe_closure_claimed"])
        self.assertTrue(all(self.payload["gatekeeper_checks"].values()))

    def test_rg_audit_and_docs(self) -> None:
        self.assertEqual(self.payload["rg_non_duplication_audit"]["tool"], "rg")
        self.assertIn("P2566/S1516", MD.read_text(encoding="utf-8"))
        self.assertIn("P2566/S1516", (ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md").read_text(encoding="utf-8"))
        self.assertIn("P2566/S1516", (ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md").read_text(encoding="utf-8"))


if __name__ == "__main__":
    unittest.main()
