from __future__ import annotations

import json
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT))

from p2501_s1451_phase_normalized_curvature_interval_newton_root_enclosure_certificate import (
    MD,
    OUT,
    build_payload,
    write_markdown,
)


class P2501PhaseNormalizedCurvatureIntervalNewtonRootEnclosureCertificateTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.payload = build_payload()
        OUT.write_text(json.dumps(cls.payload, indent=2, sort_keys=True), encoding="utf-8")
        write_markdown(cls.payload)
        cls.theorem = cls.payload["phase_normalized_curvature_interval_newton_root_enclosure_certificate"]["theorem_export"]
        cls.cert = cls.theorem["interval_newton_root_enclosure_certificate"]
        cls.step = cls.cert["interval_newton_step"]

    def test_identity_and_status(self) -> None:
        self.assertEqual(self.payload["packet_id"], "P2501")
        self.assertEqual(self.payload["stage_id"], "S1451")
        self.assertIn("INTERVAL_NEWTON_ROOT_ENCLOSURE", self.payload["status"])
        self.assertIn("NO_FORMAL_DIRECTED_BACKEND", self.payload["status"])
        self.assertIn("NO_GLOBAL_ANALYTIC_UNIQUENESS", self.payload["status"])
        self.assertIn("NO_SOURCE_EXPORT", self.payload["status"])
        self.assertIn("NO_BRIDGE_ATOM", self.payload["status"])
        self.assertIn("NO_QW2191", self.payload["status"])
        self.assertIn("NO_ROLE_TRANSFER", self.payload["status"])
        self.assertIn("NO_TOE", self.payload["status"])

    def test_interval_newton_contraction(self) -> None:
        self.assertEqual(self.theorem["starting_window_d"], ["0.34961674", "0.34961675"])
        start_left, start_right = [float(value) for value in self.theorem["starting_window_d"]]
        contracted_left, contracted_right = [float(value) for value in self.theorem["contracted_root_enclosure_d"]]
        self.assertGreater(contracted_left, start_left)
        self.assertLess(contracted_right, start_right)
        self.assertLess(float(self.theorem["contracted_root_enclosure_width"]), float(self.theorem["starting_window_width"]))
        self.assertGreater(float(self.theorem["contraction_factor"]), 1.0e7)
        self.assertTrue(self.theorem["newton_image_subset_of_starting_window"])
        self.assertTrue(self.theorem["newton_image_strict_subset_of_starting_window"])

    def test_derivative_and_endpoint_signs(self) -> None:
        self.assertEqual(self.theorem["third_derivative_interval_sign"], -1)
        self.assertTrue(self.theorem["third_derivative_strictly_negative"])
        self.assertEqual(self.step["left_endpoint_curvature_interval"]["curvature_interval_sign"], 1)
        self.assertEqual(self.step["right_endpoint_curvature_interval"]["curvature_interval_sign"], -1)
        self.assertTrue(self.theorem["contracted_endpoint_intervals_have_opposite_strict_signs"])
        self.assertTrue(self.theorem["contractive_interval_newton_enclosure_certified"])
        self.assertTrue(self.theorem["unique_inflection_root_in_contracted_enclosure_on_audited_branch"])

    def test_inherited_guards_and_gatekeepers(self) -> None:
        self.assertTrue(self.theorem["p2499_local_uniqueness_inherited"])
        self.assertTrue(self.theorem["p2499_outside_refined_window_zero_exclusion_inherited"])
        self.assertTrue(self.theorem["p2500_symbolic_formula_provenance_inherited"])
        self.assertTrue(self.theorem["single_inflection_root_on_audited_domain_inherited_and_narrowed"])
        self.assertTrue(all(self.payload["gatekeeper_checks"].values()))

    def test_negative_controls(self) -> None:
        self.assertFalse(self.theorem["formal_directed_rounding_backend_exported"])
        self.assertFalse(self.theorem["global_analytic_inflection_uniqueness_theorem_exported"])
        self.assertFalse(self.theorem["curvature_dynamic_source_exported"])
        self.assertFalse(self.theorem["legacy_to_strict_bridge_atom_exported"])
        self.assertFalse(self.theorem["strict_compression_dynamic_source_exported"])
        self.assertFalse(self.theorem["selector_source_theorem_exported"])
        self.assertFalse(self.theorem["qw2191_discharged_by_this_certificate"])
        self.assertFalse(self.theorem["role_transfer_licensed_by_this_certificate"])
        self.assertFalse(self.theorem["toe_closure_exported"])

    def test_rg_audit_and_documentation_markers(self) -> None:
        audit = self.payload["rg_non_duplication_audit"]
        self.assertEqual(audit["tool"], "rg")
        self.assertIn("new_packet", audit["patterns"])
        self.assertGreaterEqual(audit["patterns"]["new_packet"]["count"], 0)
        self.assertIn("P2501/S1451", MD.read_text(encoding="utf-8"))
        self.assertIn("P2501/S1451", (ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md").read_text(encoding="utf-8"))
        self.assertIn("P2501/S1451", (ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md").read_text(encoding="utf-8"))


if __name__ == "__main__":
    unittest.main()
