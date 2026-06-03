from __future__ import annotations

import json
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT))

from p2511_s1461_strict_damping_rg_natural_spline_collapse_certificate import (
    MD,
    OUT,
    build_payload,
    write_markdown,
)


class P2511StrictDampingRGNaturalSplineCollapseCertificateTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.payload = build_payload()
        OUT.write_text(json.dumps(cls.payload, indent=2, sort_keys=True), encoding="utf-8")
        write_markdown(cls.payload)
        cls.theorem = cls.payload["strict_damping_rg_natural_spline_collapse_certificate"]["theorem_export"]
        cls.cert = cls.theorem["strict_damping_rg_natural_spline_collapse_certificate"]

    def test_identity_and_status(self) -> None:
        self.assertEqual(self.payload["packet_id"], "P2511")
        self.assertEqual(self.payload["stage_id"], "S1461")
        self.assertIn("NATURAL_SPLINE_COLLAPSE", self.payload["status"])
        self.assertIn("NO_SOURCE_EXPORT", self.payload["status"])
        self.assertIn("NO_BRIDGE_THEOREM", self.payload["status"])
        self.assertIn("NO_ROLE_TRANSFER", self.payload["status"])
        self.assertIn("NO_QW2191", self.payload["status"])
        self.assertIn("NO_TOE", self.payload["status"])

    def test_symbolic_and_spline_collapse(self) -> None:
        self.assertEqual(self.theorem["frontier_atom_under_attack"], "strict_damping_beta_eta_source")
        self.assertTrue(self.theorem["p2510_kkt_stationarity_inherited"])
        self.assertTrue(self.theorem["symbolic_slope_residuals_zero"])
        self.assertTrue(self.theorem["tridiagonal_leading_minors_positive"])
        self.assertTrue(self.theorem["tridiagonal_pivots_positive"])
        self.assertLess(float(self.theorem["max_abs_rhs"]), 1e-90)
        self.assertLess(float(self.theorem["max_abs_second_derivative_knot"]), 1e-90)
        self.assertLess(float(self.theorem["max_abs_quadratic_coefficient"]), 1e-90)
        self.assertLess(float(self.theorem["max_abs_cubic_coefficient"]), 1e-90)
        self.assertLess(float(self.theorem["closed_form_roughness_energy"]), 1e-180)
        self.assertTrue(self.theorem["natural_spline_collapse_theorem_for_postulated_functional"])

    def test_conditional_status_and_negative_controls(self) -> None:
        self.assertTrue(self.theorem["roughness_action_still_postulated_not_derived"])
        self.assertFalse(self.theorem["strict_damping_beta_eta_source_exported"])
        self.assertFalse(self.theorem["strict_dynamical_source_for_A_P_D_exported"])
        self.assertFalse(self.theorem["strict_phase_frequency_source_exported"])
        self.assertFalse(self.theorem["bridge_theorem_exported"])
        self.assertFalse(self.theorem["role_transfer_theorem_exported"])
        self.assertFalse(self.theorem["selector_closure_exported"])
        self.assertFalse(self.theorem["qw2191_discharged_by_this_certificate"])
        self.assertFalse(self.theorem["toe_closure_claimed"])
        self.assertTrue(all(self.payload["gatekeeper_checks"].values()))

    def test_rg_audit_and_documentation_markers(self) -> None:
        audit = self.payload["rg_non_duplication_audit"]
        self.assertEqual(audit["tool"], "rg")
        self.assertIn("new_packet", audit["patterns"])
        self.assertGreaterEqual(audit["patterns"]["new_packet"]["count"], 0)
        self.assertIn("P2511/S1461", MD.read_text(encoding="utf-8"))
        self.assertIn("P2511/S1461", (ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md").read_text(encoding="utf-8"))
        self.assertIn("P2511/S1461", (ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md").read_text(encoding="utf-8"))


if __name__ == "__main__":
    unittest.main()
