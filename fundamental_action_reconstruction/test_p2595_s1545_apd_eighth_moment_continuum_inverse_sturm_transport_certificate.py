from __future__ import annotations

import json
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT))

from p2595_s1545_apd_eighth_moment_continuum_inverse_sturm_transport_certificate import MD, OUT, append_doc_sections, build_payload, write_markdown


class P2595APDEighthMomentContinuumInverseSturmTransportCertificateTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.payload = build_payload()
        OUT.write_text(json.dumps(cls.payload, indent=2, sort_keys=True), encoding="utf-8")
        write_markdown(cls.payload)
        append_doc_sections()
        cls.theorem = cls.payload["apd_eighth_moment_continuum_inverse_sturm_transport_certificate"]["theorem_export"]
        cls.certificate = cls.theorem["apd_eighth_moment_continuum_inverse_sturm_transport_certificate"]

    def test_identity_and_inheritance(self) -> None:
        self.assertEqual(self.payload["packet_id"], "P2595")
        self.assertEqual(self.payload["stage_id"], "S1545")
        self.assertIn("EIGHTH_MOMENT_CONTINUUM_INVERSE_STURM_TRANSPORT", self.payload["status"])
        self.assertEqual(self.theorem["frontier_atom_under_attack"], "strict_dynamical_source_for_A_P_D")
        self.assertTrue(self.theorem["p2591_sturm_interval_inherited"])
        self.assertTrue(self.theorem["p2594_conditional_inverse_inherited"])

    def test_continuum_inverse_transport(self) -> None:
        self.assertEqual(self.certificate["product_parameter_interval"], [300, 576])
        self.assertEqual(self.certificate["internal_eighth_interval_ascending_exact"], ["72354", "73458"])
        self.assertEqual(self.certificate["central_eighth_interval_ascending_exact"], ["42564651105/32768", "42637002849/32768"])
        self.assertEqual(self.certificate["internal_inverse_formula"], "37329/2 - p4/4")
        self.assertEqual(self.certificate["central_inverse_formula"], "42715646049/262144 - m8/8")
        self.assertTrue(self.theorem["continuum_internal_eighth_interval_maps_to_valid_product_interval"])
        self.assertTrue(self.theorem["continuum_central_eighth_interval_maps_to_valid_product_interval"])
        self.assertTrue(self.theorem["p2591_sturm_validity_transports_to_entire_eighth_shell_interval"])
        self.assertEqual(self.certificate["internal_affine_slope"], -4.0)
        self.assertEqual(self.certificate["central_affine_slope"], -8.0)
        self.assertEqual(self.certificate["internal_inverse_slope"], -0.25)
        self.assertEqual(self.certificate["central_inverse_slope"], -0.125)
        for snapshot in [*self.certificate["endpoint_support_snapshots"], self.certificate["midpoint_support_snapshot"]]:
            self.assertEqual(snapshot["point_count"], 10)
            self.assertEqual(snapshot["min_point"], 0.25)
            self.assertEqual(snapshot["max_point"], 10.75)

    def test_negative_controls_and_docs(self) -> None:
        self.assertFalse(self.theorem["apd_eighth_moment_continuum_source_exported"])
        self.assertFalse(self.theorem["apd_continuum_inverse_selector_source_exported"])
        self.assertFalse(self.theorem["apd_inverse_sturm_transport_source_exported"])
        self.assertFalse(self.theorem["apd_continuum_support_reconstruction_source_exported"])
        self.assertFalse(self.theorem["strict_dynamical_source_for_A_P_D_exported"])
        self.assertFalse(self.theorem["bridge_theorem_exported"])
        self.assertFalse(self.theorem["role_transfer_theorem_exported"])
        self.assertFalse(self.theorem["qw2191_discharged_by_this_certificate"])
        self.assertFalse(self.theorem["toe_closure_claimed"])
        self.assertTrue(all(self.payload["gatekeeper_checks"].values()))
        self.assertEqual(self.payload["rg_non_duplication_audit"]["tool"], "rg")
        self.assertIn("P2595/S1545", MD.read_text(encoding="utf-8"))
        self.assertIn("P2595/S1545", (ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md").read_text(encoding="utf-8"))
        self.assertIn("P2595/S1545", (ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md").read_text(encoding="utf-8"))


if __name__ == "__main__":
    unittest.main()
