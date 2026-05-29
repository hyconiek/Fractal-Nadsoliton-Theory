from __future__ import annotations

import hashlib
import json
import subprocess
import sys
import unittest
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parent
SCRIPT = ROOT / "p2351_s1301_strict_track_b_quartic_support_global_convexity_jacobian_identity_witness.py"
OUT = ROOT / "generated" / "p2351_s1301_strict_track_b_quartic_support_global_convexity_jacobian_identity_witness.json"
MD = ROOT / "generated" / "p2351_s1301_strict_track_b_quartic_support_global_convexity_jacobian_identity_witness.md"


def sha256_json(obj: Any) -> str:
    return hashlib.sha256(
        json.dumps(obj, sort_keys=True, ensure_ascii=False, separators=(",", ":")).encode("utf-8")
    ).hexdigest()


class P2351QuarticSupportGlobalConvexityJacobianTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        subprocess.run([sys.executable, str(SCRIPT)], check=True, cwd=ROOT.parent)
        cls.payload = json.loads(OUT.read_text(encoding="utf-8"))
        cls.probe = cls.payload["strict_track_b_quartic_support_global_convexity_jacobian_identity_witness_probe"]
        cls.witness = cls.probe["track_B_quartic_support_global_convexity_jacobian_identity_witness"]
        cls.convexity = cls.witness["convexity_certificate"]
        cls.jacobian = cls.witness["jacobian_identity_certificate"]

    def test_packet_identity_and_outputs(self) -> None:
        self.assertEqual(self.payload["schema_version"], "p2351_s1301_v1")
        self.assertEqual(self.payload["packet_id"], "P2351")
        self.assertEqual(self.payload["stage_id"], "S1301")
        self.assertTrue(OUT.exists())
        self.assertTrue(MD.exists())

    def test_global_convexity_certificate(self) -> None:
        self.assertEqual(self.witness["ledger_coefficient_b_GB"], "13152087*log(2)/(320000000*pi**2)")
        self.assertEqual(self.witness["per_euler_pairing"], "13152087*log(2)/10000000")
        radii = self.convexity["radii_eigenvalues"]
        self.assertEqual(radii["r_perp_multiplicity_2"], "-(3*y**2 - 10)/10")
        self.assertEqual(radii["r_gradient_multiplicity_1"], "-(15*y**2 - 12*y - 10)/10")
        self.assertEqual(self.convexity["min_radius_lower_bound"], "7/10")
        self.assertEqual(self.convexity["positive_certificates_on_0_1"]["r_perp_minus_7_over_10"], "-3*(y - 1)*(y + 1)/10")
        self.assertEqual(
            self.convexity["positive_certificates_on_0_1"]["r_gradient_minus_7_over_10"],
            "-3*(y - 1)*(5*y + 1)/10",
        )
        self.assertTrue(self.convexity["global_strict_convexity_for_eps_1_over_10"])
        self.assertEqual(self.convexity["determinant_lower_bound_from_radii"], "343/1000")

    def test_jacobian_identity_integrated_replay_and_gatekeepers(self) -> None:
        self.assertEqual(self.jacobian["density_times_jacobian_minus_16_residual"], "0")
        self.assertEqual(self.jacobian["abstract_q_residual"], "0")
        self.assertEqual(self.jacobian["integral_sigma3_dA_via_gauss_map_degree_1"], "2*pi**2")
        self.assertEqual(self.jacobian["integrated_boundary_functional_16_integral_sigma3"], "32*pi**2")
        self.assertEqual(self.jacobian["integrated_boundary_residual"], "0")
        self.assertEqual(self.jacobian["normalized_track_B_pairing"], "13152087*log(2)/10000000")
        self.assertEqual(self.jacobian["normalized_pairing_residual"], "0")
        self.assertEqual(self.jacobian["p2350_integrated_replay_residual"], "0")
        self.assertEqual(self.jacobian["p2350_pairing_replay_residual"], "0")
        for key, value in self.payload["gatekeeper_checks"].items():
            self.assertTrue(value, key)
        theorem = self.probe["theorem_export"]
        self.assertEqual(self.probe["theorem_fingerprint_sha256"], sha256_json(theorem))
        self.assertIn("general support-function theorem for arbitrary h on S^3", theorem["not_licensed"])
        self.assertIn("selector premise or QW-2191 selector discharge", theorem["not_licensed"])
        self.assertIn("ToE closure", theorem["not_licensed"])


if __name__ == "__main__":
    unittest.main()
