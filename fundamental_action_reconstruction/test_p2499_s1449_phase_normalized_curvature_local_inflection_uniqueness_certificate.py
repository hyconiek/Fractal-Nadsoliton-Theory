from __future__ import annotations

import json
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT))

from p2499_s1449_phase_normalized_curvature_local_inflection_uniqueness_certificate import (
    MD,
    OUT,
    build_payload,
    write_markdown,
)


class P2499PhaseNormalizedCurvatureLocalInflectionUniquenessCertificateTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.payload = build_payload()
        OUT.write_text(json.dumps(cls.payload, indent=2, sort_keys=True), encoding="utf-8")
        write_markdown(cls.payload)
        cls.theorem = cls.payload["phase_normalized_curvature_local_inflection_uniqueness_certificate"]["theorem_export"]
        cls.cert = cls.theorem["local_inflection_uniqueness_certificate"]
        cls.third = cls.cert["third_derivative_interval_certificate"]

    def test_identity_and_status(self) -> None:
        self.assertEqual(self.payload["packet_id"], "P2499")
        self.assertEqual(self.payload["stage_id"], "S1449")
        self.assertIn("LOCAL_INFLECTION_UNIQUENESS", self.payload["status"])
        self.assertIn("NO_FORMAL_DIRECTED_BACKEND", self.payload["status"])
        self.assertIn("NO_GLOBAL_ANALYTIC_UNIQUENESS", self.payload["status"])
        self.assertIn("NO_SOURCE_EXPORT", self.payload["status"])
        self.assertIn("NO_BRIDGE_ATOM", self.payload["status"])
        self.assertIn("NO_QW2191", self.payload["status"])
        self.assertIn("NO_ROLE_TRANSFER", self.payload["status"])
        self.assertIn("NO_TOE", self.payload["status"])

    def test_refined_window_and_endpoint_signs(self) -> None:
        self.assertEqual(self.theorem["refined_inflection_window_d"], ["0.34961674", "0.34961675"])
        self.assertLessEqual(float(self.theorem["refined_inflection_window_width"]), 1e-8)
        self.assertTrue(self.theorem["endpoint_intervals_have_opposite_strict_signs"])
        self.assertEqual(self.cert["left_endpoint_curvature_interval"]["curvature_interval_sign"], 1)
        self.assertEqual(self.cert["right_endpoint_curvature_interval"]["curvature_interval_sign"], -1)

    def test_third_derivative_interval_monotonicity(self) -> None:
        bounds = self.theorem["x_third_interval_over_refined_window"]
        self.assertEqual(self.theorem["x_third_interval_sign"], -1)
        self.assertLess(float(bounds[0]), 0.0)
        self.assertLess(float(bounds[1]), 0.0)
        self.assertTrue(self.third["x_third_strictly_negative_on_refined_window"])
        self.assertTrue(self.theorem["x_second_strictly_decreasing_on_refined_window"])
        self.assertTrue(self.theorem["finite_local_unique_curvature_zero_in_refined_window"])
        self.assertTrue(self.theorem["single_local_inflection_plus_outside_exclusion_on_audited_domain"])

    def test_inherited_p2498_guard_and_gatekeepers(self) -> None:
        self.assertTrue(self.theorem["p2498_outside_refined_window_zero_exclusion_inherited"])
        self.assertTrue(all(self.payload["gatekeeper_checks"].values()))

    def test_negative_controls(self) -> None:
        self.assertFalse(self.theorem["formal_directed_rounding_backend_exported"])
        self.assertFalse(self.theorem["global_analytic_inflection_uniqueness_theorem_exported"])
        self.assertFalse(self.theorem["analytic_third_derivative_symbolic_proof_exported"])
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
        self.assertIn("P2499/S1449", MD.read_text(encoding="utf-8"))
        self.assertIn("P2499/S1449", (ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md").read_text(encoding="utf-8"))
        self.assertIn("P2499/S1449", (ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md").read_text(encoding="utf-8"))


if __name__ == "__main__":
    unittest.main()
