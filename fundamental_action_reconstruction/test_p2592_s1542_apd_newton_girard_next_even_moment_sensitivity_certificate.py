from __future__ import annotations

import json
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT))

from p2590_s1540_apd_finite_even_moment_shell_interval_nonuniqueness_audit import PRODUCT_PARAMETER_GRID
from p2592_s1542_apd_newton_girard_next_even_moment_sensitivity_certificate import MD, OUT, append_doc_sections, build_payload, write_markdown


class P2592APDNewtonGirardNextEvenMomentSensitivityCertificateTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.payload = build_payload()
        OUT.write_text(json.dumps(cls.payload, indent=2, sort_keys=True), encoding="utf-8")
        write_markdown(cls.payload)
        append_doc_sections()
        cls.theorem = cls.payload["apd_newton_girard_next_even_moment_sensitivity_certificate"]["theorem_export"]
        cls.certificate = cls.theorem["apd_newton_girard_next_even_moment_sensitivity_certificate"]

    def test_identity_and_inheritance(self) -> None:
        self.assertEqual(self.payload["packet_id"], "P2592")
        self.assertEqual(self.payload["stage_id"], "S1542")
        self.assertIn("NEWTON_GIRARD_NEXT_EVEN_MOMENT", self.payload["status"])
        self.assertEqual(self.theorem["frontier_atom_under_attack"], "strict_dynamical_source_for_A_P_D")
        self.assertTrue(self.theorem["p2591_sturm_interval_inherited"])

    def test_newton_girard_next_shell_sensitivity(self) -> None:
        self.assertEqual(self.certificate["newton_girard_formula_internal_eighth_shell_p4"], "74658.0 - 4*e4")
        self.assertEqual(self.certificate["internal_eighth_shell_slope_by_product_parameter"], -4.0)
        self.assertEqual(self.certificate["central_eighth_moment_slope_by_product_parameter"], -8.0)
        self.assertEqual(self.certificate["grid_internal_eighth_shell_distinct_count"], len(PRODUCT_PARAMETER_GRID))
        self.assertEqual(self.certificate["grid_central_eighth_moment_distinct_count"], len(PRODUCT_PARAMETER_GRID))
        self.assertLessEqual(self.certificate["max_formula_numeric_abs_error"], 1.0e-6)
        self.assertTrue(self.certificate["lower_shells_constant_but_next_shell_varies"])
        for row in self.certificate["selector_recovery_from_next_moment"]:
            self.assertAlmostEqual(row["product_parameter"], row["recovered_product_from_internal_eighth_shell"], places=8)
            self.assertAlmostEqual(row["product_parameter"], row["recovered_product_from_central_eighth_moment"], places=8)
        for snapshot in self.certificate["grid_snapshots"]:
            self.assertTrue(snapshot["has_expected_internal_second_fourth_sixth_shell"])
            self.assertEqual(snapshot["vandermonde_rank"], 10)

    def test_negative_controls_and_docs(self) -> None:
        self.assertFalse(self.theorem["apd_eighth_moment_shell_source_exported"])
        self.assertFalse(self.theorem["apd_next_even_moment_selector_source_exported"])
        self.assertFalse(self.theorem["apd_newton_girard_selector_source_exported"])
        self.assertFalse(self.theorem["apd_product_parameter_source_exported"])
        self.assertFalse(self.theorem["strict_dynamical_source_for_A_P_D_exported"])
        self.assertFalse(self.theorem["bridge_theorem_exported"])
        self.assertFalse(self.theorem["role_transfer_theorem_exported"])
        self.assertFalse(self.theorem["qw2191_discharged_by_this_certificate"])
        self.assertFalse(self.theorem["toe_closure_claimed"])
        self.assertTrue(all(self.payload["gatekeeper_checks"].values()))
        self.assertEqual(self.payload["rg_non_duplication_audit"]["tool"], "rg")
        self.assertIn("P2592/S1542", MD.read_text(encoding="utf-8"))
        self.assertIn("P2592/S1542", (ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md").read_text(encoding="utf-8"))
        self.assertIn("P2592/S1542", (ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md").read_text(encoding="utf-8"))


if __name__ == "__main__":
    unittest.main()
