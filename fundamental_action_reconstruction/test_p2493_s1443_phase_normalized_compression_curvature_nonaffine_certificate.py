from __future__ import annotations

import json
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT))

from p2493_s1443_phase_normalized_compression_curvature_nonaffine_certificate import (
    MD,
    OUT,
    build_payload,
    write_markdown,
)


class P2493PhaseNormalizedCompressionCurvatureNonaffineCertificateTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.payload = build_payload()
        OUT.write_text(json.dumps(cls.payload, indent=2, sort_keys=True), encoding="utf-8")
        write_markdown(cls.payload)
        cls.theorem = cls.payload["phase_normalized_compression_curvature_nonaffine_certificate"]["theorem_export"]
        cls.cert = cls.theorem["phase_normalized_compression_curvature_certificate"]

    def test_identity_and_status(self) -> None:
        self.assertEqual(self.payload["packet_id"], "P2493")
        self.assertEqual(self.payload["stage_id"], "S1443")
        self.assertIn("COMPRESSION_CURVATURE_NONAFFINE", self.payload["status"])
        self.assertIn("NO_SOURCE_EXPORT", self.payload["status"])
        self.assertIn("NO_BRIDGE_ATOM", self.payload["status"])
        self.assertIn("NO_QW2191", self.payload["status"])
        self.assertIn("NO_ROLE_TRANSFER", self.payload["status"])
        self.assertIn("NO_TOE", self.payload["status"])

    def test_curvature_counts_and_identity(self) -> None:
        self.assertEqual(self.theorem["sample_count"], 10)
        self.assertEqual(self.theorem["nonzero_curvature_sample_count"], 10)
        self.assertLess(self.theorem["max_second_derivative_identity_residual_abs"], 1e-12)
        self.assertTrue(self.theorem["has_positive_curvature_near_origin"])
        self.assertTrue(self.theorem["has_negative_curvature_on_tail"])
        self.assertTrue(self.theorem["affine_bridge_ruled_out_by_curvature"])
        self.assertEqual(self.cert["sample_curvature_sign_sequence"][:4], [1, 1, 1, 1])
        self.assertEqual(self.cert["sample_curvature_sign_sequence"][4:], [-1, -1, -1, -1, -1, -1])

    def test_inflection_and_z12_discrete_curvature(self) -> None:
        inflection = self.cert["inflection_certificate"]
        self.assertEqual(inflection["bracket"], [0.1, 0.5])
        self.assertTrue(inflection["sign_change_certified"])
        self.assertGreater(inflection["root_estimate"], 0.34)
        self.assertLess(inflection["root_estimate"], 0.36)
        self.assertLess(inflection["root_curvature_abs"], 1e-12)
        z12 = self.cert["discrete_z12_curvature"]
        self.assertTrue(z12["all_tail_second_differences_negative_from_d1"])
        self.assertEqual(len(z12["x_values"]), 12)
        self.assertEqual(len(z12["second_differences_d1_to_d10"]), 10)
        self.assertAlmostEqual(z12["max_abs_second_difference"], 0.45966570666655304, places=12)
        self.assertTrue(all(value < 0.0 for value in z12["second_differences_d1_to_d10"]))

    def test_upstream_obligations_and_negative_controls(self) -> None:
        self.assertTrue(self.theorem["p2414_damping_nonabsorption_inherited"])
        self.assertTrue(self.theorem["p2492_shared_core_still_requires_damping_atoms"])
        self.assertTrue(all(self.payload["gatekeeper_checks"].values()))
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
        self.assertIn("P2493/S1443", MD.read_text(encoding="utf-8"))
        self.assertIn("P2493/S1443", (ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md").read_text(encoding="utf-8"))
        self.assertIn("P2493/S1443", (ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md").read_text(encoding="utf-8"))


if __name__ == "__main__":
    unittest.main()
