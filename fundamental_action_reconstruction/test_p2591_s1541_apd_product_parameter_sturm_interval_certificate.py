from __future__ import annotations

import json
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT))

from p2591_s1541_apd_product_parameter_sturm_interval_certificate import MD, OUT, append_doc_sections, build_payload, write_markdown


class P2591APDProductParameterSturmIntervalCertificateTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.payload = build_payload()
        OUT.write_text(json.dumps(cls.payload, indent=2, sort_keys=True), encoding="utf-8")
        write_markdown(cls.payload)
        append_doc_sections()
        cls.theorem = cls.payload["apd_product_parameter_sturm_interval_certificate"]["theorem_export"]
        cls.certificate = cls.theorem["apd_product_parameter_sturm_interval_certificate"]

    def test_identity_and_inheritance(self) -> None:
        self.assertEqual(self.payload["packet_id"], "P2591")
        self.assertEqual(self.payload["stage_id"], "S1541")
        self.assertIn("STURM_INTERVAL", self.payload["status"])
        self.assertEqual(self.theorem["frontier_atom_under_attack"], "strict_dynamical_source_for_A_P_D")
        self.assertTrue(self.theorem["p2590_finite_even_moment_shell_interval_inherited"])

    def test_sturm_interval_certificate(self) -> None:
        self.assertEqual(self.certificate["certified_product_parameter_interval"], [300, 576])
        self.assertIn("16*e4**3", self.certificate["discriminant"])
        self.assertEqual(self.certificate["discriminant_sturm_count_on_interval"]["root_count_open_interval"], 0)
        self.assertTrue(self.certificate["discriminant_positive_at_left_endpoint"])
        self.assertTrue(self.certificate["discriminant_positive_at_right_endpoint"])
        self.assertEqual(self.certificate["positive_root_count_anchor_e4_300"]["root_count_open_interval"], 4)
        self.assertEqual(self.certificate["negative_root_count_anchor_e4_300"]["root_count_open_interval"], 0)
        self.assertTrue(self.certificate["positive_real_root_count_is_constant_on_interval"])
        self.assertEqual(self.certificate["vieta_power_sums_fixed_on_interval"]["p1_internal_second_shell"], 30.0)
        self.assertEqual(self.certificate["vieta_power_sums_fixed_on_interval"]["p2_internal_fourth_shell"], 354.0)
        self.assertEqual(self.certificate["vieta_power_sums_fixed_on_interval"]["p3_internal_sixth_shell"], 4890.0)
        for snapshot in self.certificate["endpoint_support_snapshots"]:
            self.assertTrue(snapshot["has_expected_internal_second_fourth_sixth_shell"])
            self.assertEqual(snapshot["vandermonde_rank"], 10)
        self.assertTrue(self.theorem["continuous_interval_of_valid_supports_certified"])
        self.assertTrue(self.theorem["finite_even_moment_shell_interval_still_does_not_select_apd_support"])

    def test_negative_controls_and_docs(self) -> None:
        self.assertFalse(self.theorem["apd_product_parameter_interval_source_exported"])
        self.assertFalse(self.theorem["apd_sturm_interval_selector_source_exported"])
        self.assertFalse(self.theorem["apd_finite_even_moment_shell_source_exported"])
        self.assertFalse(self.theorem["apd_support_selection_source_exported"])
        self.assertFalse(self.theorem["strict_dynamical_source_for_A_P_D_exported"])
        self.assertFalse(self.theorem["bridge_theorem_exported"])
        self.assertFalse(self.theorem["role_transfer_theorem_exported"])
        self.assertFalse(self.theorem["qw2191_discharged_by_this_certificate"])
        self.assertFalse(self.theorem["toe_closure_claimed"])
        self.assertTrue(all(self.payload["gatekeeper_checks"].values()))
        self.assertEqual(self.payload["rg_non_duplication_audit"]["tool"], "rg")
        self.assertIn("P2591/S1541", MD.read_text(encoding="utf-8"))
        self.assertIn("P2591/S1541", (ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md").read_text(encoding="utf-8"))
        self.assertIn("P2591/S1541", (ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md").read_text(encoding="utf-8"))


if __name__ == "__main__":
    unittest.main()
