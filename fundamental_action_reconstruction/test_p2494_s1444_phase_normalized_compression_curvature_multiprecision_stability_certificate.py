from __future__ import annotations

import json
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT))

from p2494_s1444_phase_normalized_compression_curvature_multiprecision_stability_certificate import (
    MD,
    OUT,
    build_payload,
    write_markdown,
)


class P2494PhaseNormalizedCompressionCurvatureMultiprecisionStabilityCertificateTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.payload = build_payload()
        OUT.write_text(json.dumps(cls.payload, indent=2, sort_keys=True), encoding="utf-8")
        write_markdown(cls.payload)
        cls.theorem = cls.payload["phase_normalized_compression_curvature_multiprecision_stability_certificate"]["theorem_export"]
        cls.cert = cls.theorem["multiprecision_stability_certificate"]

    def test_identity_and_status(self) -> None:
        self.assertEqual(self.payload["packet_id"], "P2494")
        self.assertEqual(self.payload["stage_id"], "S1444")
        self.assertIn("MULTIPRECISION_STABILITY", self.payload["status"])
        self.assertIn("NO_DIRECTED_ROUNDING", self.payload["status"])
        self.assertIn("NO_SOURCE_EXPORT", self.payload["status"])
        self.assertIn("NO_BRIDGE_ATOM", self.payload["status"])
        self.assertIn("NO_QW2191", self.payload["status"])
        self.assertIn("NO_ROLE_TRANSFER", self.payload["status"])
        self.assertIn("NO_TOE", self.payload["status"])

    def test_precision_levels_and_stability_flags(self) -> None:
        self.assertEqual(self.theorem["dps_levels"], [50, 80, 120])
        self.assertEqual([row["dps"] for row in self.cert["precision_rows"]], [50, 80, 120])
        self.assertTrue(self.theorem["sample_sign_sequence_stable_across_precisions"])
        self.assertTrue(self.theorem["z12_sign_sequence_stable_across_precisions"])
        self.assertTrue(self.theorem["nonzero_count_stable"])
        self.assertTrue(self.theorem["inflection_sign_change_stable"])
        self.assertTrue(self.theorem["p2493_affine_bridge_ruled_out_inherited"])

    def test_highest_precision_margins_and_signs(self) -> None:
        highest = self.cert["precision_rows"][-1]
        self.assertEqual(highest["sample_curvature_sign_sequence"][:4], [1, 1, 1, 1])
        self.assertEqual(highest["sample_curvature_sign_sequence"][4:], [-1, -1, -1, -1, -1, -1])
        self.assertEqual(highest["nonzero_curvature_sample_count"], 10)
        self.assertGreater(float(self.theorem["highest_precision_min_abs_sample_curvature"]), 6e-4)
        self.assertGreater(float(self.theorem["highest_precision_min_abs_z12_second_difference"]), 9e-4)
        self.assertTrue(highest["inflection_sign_change_certified"])
        self.assertTrue(highest["z12_summary"]["all_tail_second_differences_negative"])
        self.assertEqual(highest["z12_summary"]["second_difference_sign_sequence"], [-1] * 10)

    def test_drift_recording_and_gatekeepers(self) -> None:
        self.assertGreaterEqual(float(self.theorem["max_x_second_drift_80_to_120"]), 0.0)
        self.assertGreaterEqual(float(self.theorem["inflection_root_drift_80_to_120"]), 0.0)
        self.assertTrue(all(self.payload["gatekeeper_checks"].values()))

    def test_negative_controls(self) -> None:
        self.assertFalse(self.theorem["directed_rounding_interval_proof_exported"])
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
        self.assertIn("P2494/S1444", MD.read_text(encoding="utf-8"))
        self.assertIn("P2494/S1444", (ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md").read_text(encoding="utf-8"))
        self.assertIn("P2494/S1444", (ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md").read_text(encoding="utf-8"))


if __name__ == "__main__":
    unittest.main()
