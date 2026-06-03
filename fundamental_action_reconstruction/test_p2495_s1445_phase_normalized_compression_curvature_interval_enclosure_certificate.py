from __future__ import annotations

import json
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT))

from p2495_s1445_phase_normalized_compression_curvature_interval_enclosure_certificate import (
    MD,
    OUT,
    build_payload,
    write_markdown,
)


class P2495PhaseNormalizedCompressionCurvatureIntervalEnclosureCertificateTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.payload = build_payload()
        OUT.write_text(json.dumps(cls.payload, indent=2, sort_keys=True), encoding="utf-8")
        write_markdown(cls.payload)
        cls.theorem = cls.payload["phase_normalized_compression_curvature_interval_enclosure_certificate"]["theorem_export"]
        cls.cert = cls.theorem["interval_enclosure_certificate"]

    def test_identity_and_status(self) -> None:
        self.assertEqual(self.payload["packet_id"], "P2495")
        self.assertEqual(self.payload["stage_id"], "S1445")
        self.assertIn("INTERVAL_ENCLOSURE", self.payload["status"])
        self.assertIn("NO_FORMAL_DIRECTED_ROUNDING", self.payload["status"])
        self.assertIn("NO_SOURCE_EXPORT", self.payload["status"])
        self.assertIn("NO_BRIDGE_ATOM", self.payload["status"])
        self.assertIn("NO_QW2191", self.payload["status"])
        self.assertIn("NO_ROLE_TRANSFER", self.payload["status"])
        self.assertIn("NO_TOE", self.payload["status"])

    def test_interval_backend_and_sample_signs(self) -> None:
        self.assertEqual(self.theorem["interval_backend"], "mpmath.iv")
        self.assertEqual(self.theorem["interval_dps"], 100)
        self.assertEqual(self.theorem["point_solve_dps"], 120)
        self.assertEqual(self.theorem["sample_interval_sign_sequence"], [1, 1, 1, 1, -1, -1, -1, -1, -1, -1])
        self.assertTrue(self.theorem["all_sample_curvature_intervals_strictly_signed"])
        self.assertTrue(self.theorem["sample_interval_signs_match_p2493_p2494"])
        self.assertEqual(len(self.cert["sample_rows"]), 10)
        self.assertTrue(all(row["curvature_interval_strictly_signed"] for row in self.cert["sample_rows"]))

    def test_interval_bounds_have_expected_signs(self) -> None:
        first = self.cert["sample_rows"][0]
        tail = self.cert["sample_rows"][-1]
        self.assertGreater(float(first["x_second_interval"][0]), 19.0)
        self.assertGreater(float(first["x_second_interval"][1]), 19.0)
        self.assertLess(float(tail["x_second_interval"][0]), -6e-4)
        self.assertLess(float(tail["x_second_interval"][1]), -6e-4)
        self.assertGreater(float(first["x_second_interval_width"]), 0.0)

    def test_z12_interval_certificate(self) -> None:
        z12 = self.cert["z12_interval_certificate"]
        self.assertTrue(self.theorem["z12_second_difference_intervals_all_negative"])
        self.assertTrue(z12["all_second_difference_intervals_negative"])
        self.assertEqual(len(z12["rows"]), 10)
        self.assertTrue(all(row["second_difference_sign"] == -1 for row in z12["rows"]))
        self.assertGreater(float(z12["minimum_negative_margin"]), 9e-4)

    def test_inherited_stability_gatekeepers_and_negative_controls(self) -> None:
        self.assertTrue(self.theorem["p2494_sample_sign_stability_inherited"])
        self.assertTrue(self.theorem["p2494_z12_sign_stability_inherited"])
        self.assertTrue(all(self.payload["gatekeeper_checks"].values()))
        self.assertFalse(self.theorem["formal_directed_rounding_backend_exported"])
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
        self.assertIn("P2495/S1445", MD.read_text(encoding="utf-8"))
        self.assertIn("P2495/S1445", (ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md").read_text(encoding="utf-8"))
        self.assertIn("P2495/S1445", (ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md").read_text(encoding="utf-8"))


if __name__ == "__main__":
    unittest.main()
