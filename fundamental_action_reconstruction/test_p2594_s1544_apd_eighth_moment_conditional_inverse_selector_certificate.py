from __future__ import annotations

import json
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT))

from p2590_s1540_apd_finite_even_moment_shell_interval_nonuniqueness_audit import PRODUCT_PARAMETER_GRID
from p2594_s1544_apd_eighth_moment_conditional_inverse_selector_certificate import MD, OUT, append_doc_sections, build_payload, write_markdown


class P2594APDEighthMomentConditionalInverseSelectorCertificateTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.payload = build_payload()
        OUT.write_text(json.dumps(cls.payload, indent=2, sort_keys=True), encoding="utf-8")
        write_markdown(cls.payload)
        append_doc_sections()
        cls.theorem = cls.payload["apd_eighth_moment_conditional_inverse_selector_certificate"]["theorem_export"]
        cls.certificate = cls.theorem["apd_eighth_moment_conditional_inverse_selector_certificate"]

    def test_identity_and_inheritance(self) -> None:
        self.assertEqual(self.payload["packet_id"], "P2594")
        self.assertEqual(self.payload["stage_id"], "S1544")
        self.assertIn("EIGHTH_MOMENT_CONDITIONAL_INVERSE", self.payload["status"])
        self.assertEqual(self.theorem["frontier_atom_under_attack"], "strict_dynamical_source_for_A_P_D")
        self.assertTrue(self.theorem["p2593_exact_replay_inherited"])

    def test_affine_inverse_and_support_reconstruction(self) -> None:
        bounds = self.certificate["exact_interval_bounds"]
        self.assertEqual(bounds["internal_inverse_formula"], "e4 = (74658 - p4) / 4")
        self.assertEqual(bounds["central_inverse_formula"], "e4 = (42715646049/32768 - M8) / 8")
        self.assertTrue(bounds["internal_map_is_affine_bijection_on_interval"])
        self.assertTrue(bounds["central_map_is_affine_bijection_on_interval"])
        self.assertEqual(self.certificate["support_reconstruction_row_count"], len(PRODUCT_PARAMETER_GRID))
        self.assertTrue(self.certificate["internal_eighth_values_are_strictly_decreasing_on_grid"])
        self.assertTrue(self.certificate["central_eighth_values_are_strictly_decreasing_on_grid"])
        self.assertTrue(self.theorem["eighth_moment_conditionally_recovers_product_parameter"])
        self.assertTrue(self.theorem["recovered_product_conditionally_reconstructs_support_roots"])
        for row in self.certificate["support_reconstruction_rows"]:
            self.assertLessEqual(row["max_abs_recovered_product_error"], 1.0e-10)
            self.assertLessEqual(row["max_abs_squared_offset_reconstruction_error"], 1.0e-9)
            self.assertLessEqual(row["max_abs_poly_root_reconstruction_error"], 1.0e-7)
            self.assertEqual(row["support_cardinality"], 10)

    def test_negative_controls_and_docs(self) -> None:
        self.assertFalse(self.theorem["apd_eighth_moment_inverse_source_exported"])
        self.assertFalse(self.theorem["apd_conditional_support_recovery_source_exported"])
        self.assertFalse(self.theorem["apd_next_shell_conditional_selector_source_exported"])
        self.assertFalse(self.theorem["apd_root_reconstruction_source_exported"])
        self.assertFalse(self.theorem["strict_dynamical_source_for_A_P_D_exported"])
        self.assertFalse(self.theorem["bridge_theorem_exported"])
        self.assertFalse(self.theorem["role_transfer_theorem_exported"])
        self.assertFalse(self.theorem["qw2191_discharged_by_this_certificate"])
        self.assertFalse(self.theorem["toe_closure_claimed"])
        self.assertTrue(all(self.payload["gatekeeper_checks"].values()))
        self.assertEqual(self.payload["rg_non_duplication_audit"]["tool"], "rg")
        self.assertIn("P2594/S1544", MD.read_text(encoding="utf-8"))
        self.assertIn("P2594/S1544", (ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md").read_text(encoding="utf-8"))
        self.assertIn("P2594/S1544", (ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md").read_text(encoding="utf-8"))


if __name__ == "__main__":
    unittest.main()
