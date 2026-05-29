from __future__ import annotations

import hashlib
import json
import subprocess
import sys
import unittest
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parent
SCRIPT = ROOT / "p2350_s1300_strict_track_b_quartic_support_convex_boundary_integrated_witness.py"
OUT = ROOT / "generated" / "p2350_s1300_strict_track_b_quartic_support_convex_boundary_integrated_witness.json"
MD = ROOT / "generated" / "p2350_s1300_strict_track_b_quartic_support_convex_boundary_integrated_witness.md"


def sha256_json(obj: Any) -> str:
    return hashlib.sha256(
        json.dumps(obj, sort_keys=True, ensure_ascii=False, separators=(",", ":")).encode("utf-8")
    ).hexdigest()


class P2350QuarticSupportConvexBoundaryWitnessTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        subprocess.run([sys.executable, str(SCRIPT)], check=True, cwd=ROOT.parent)
        cls.payload = json.loads(OUT.read_text(encoding="utf-8"))
        cls.probe = cls.payload["strict_track_b_quartic_support_convex_boundary_integrated_witness_probe"]
        cls.witness = cls.probe["track_B_quartic_support_convex_boundary_integrated_witness"]
        cls.support = cls.witness["support_witness"]

    def test_packet_identity_and_outputs(self) -> None:
        self.assertEqual(self.payload["schema_version"], "p2350_s1300_v1")
        self.assertEqual(self.payload["packet_id"], "P2350")
        self.assertEqual(self.payload["stage_id"], "S1300")
        self.assertTrue(OUT.exists())
        self.assertTrue(MD.exists())

    def test_support_samples_and_integrated_values(self) -> None:
        self.assertEqual(self.witness["ledger_coefficient_b_GB"], "13152087*log(2)/(320000000*pi**2)")
        self.assertEqual(self.witness["per_euler_pairing"], "13152087*log(2)/10000000")
        self.assertEqual(self.support["epsilon"], "1/10")
        rows = self.support["sample_rows"]
        self.assertEqual(len(rows), 3)
        self.assertEqual(rows[0]["principal_radii"], ["7/10", "7/10", "7/10"])
        self.assertEqual(rows[0]["shape_sigma3_product_curvatures"], "1000/343")
        self.assertEqual(rows[1]["principal_radii"], ["1", "1", "1"])
        self.assertEqual(rows[2]["principal_radii"], ["193/160", "157/160", "157/160"])
        self.assertTrue(self.support["all_recorded_sample_radii_positive"])
        self.assertEqual(self.support["nonconstant_sigma3_witness_pole_minus_equator"], "657/343")
        self.assertEqual(self.support["nonconstant_density_witness_pole_minus_equator"], "10512/343")
        self.assertEqual(self.support["pointwise_density_jacobian_residual"], "0")
        self.assertEqual(self.support["integral_sigma3_via_gauss_map"], "2*pi**2")
        self.assertEqual(self.support["integrated_boundary_functional_16_integral_sigma3"], "32*pi**2")
        self.assertEqual(self.support["integrated_boundary_residual"], "0")
        self.assertEqual(self.support["normalized_track_B_pairing"], "13152087*log(2)/10000000")
        self.assertEqual(self.support["normalized_pairing_residual"], "0")

    def test_dependencies_gatekeepers_and_fingerprint(self) -> None:
        deps = self.probe["current_export_dependencies"]
        self.assertTrue(deps["p2345_boundary_value_replayed"])
        self.assertTrue(deps["p2348_flat_integral_replay_replayed"])
        self.assertTrue(deps["p2349_integrated_value_replayed"])
        for key, value in self.payload["gatekeeper_checks"].items():
            self.assertTrue(value, key)
        theorem = self.probe["theorem_export"]
        self.assertEqual(self.probe["theorem_fingerprint_sha256"], sha256_json(theorem))
        self.assertIn("integrated arbitrary-boundary theorem", theorem["not_licensed"])
        self.assertIn("selector premise or QW-2191 selector discharge", theorem["not_licensed"])
        self.assertIn("ToE closure", theorem["not_licensed"])


if __name__ == "__main__":
    unittest.main()
