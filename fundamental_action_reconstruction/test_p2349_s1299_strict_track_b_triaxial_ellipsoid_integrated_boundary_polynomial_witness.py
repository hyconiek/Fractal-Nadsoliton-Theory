from __future__ import annotations

import hashlib
import json
import subprocess
import sys
import unittest
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parent
SCRIPT = ROOT / "p2349_s1299_strict_track_b_triaxial_ellipsoid_integrated_boundary_polynomial_witness.py"
OUT = ROOT / "generated" / "p2349_s1299_strict_track_b_triaxial_ellipsoid_integrated_boundary_polynomial_witness.json"
MD = ROOT / "generated" / "p2349_s1299_strict_track_b_triaxial_ellipsoid_integrated_boundary_polynomial_witness.md"


def sha256_json(obj: Any) -> str:
    return hashlib.sha256(
        json.dumps(obj, sort_keys=True, ensure_ascii=False, separators=(",", ":")).encode("utf-8")
    ).hexdigest()


class P2349TriaxialEllipsoidIntegratedWitnessTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        subprocess.run([sys.executable, str(SCRIPT)], check=True, cwd=ROOT.parent)
        cls.payload = json.loads(OUT.read_text(encoding="utf-8"))
        cls.probe = cls.payload["strict_track_b_triaxial_ellipsoid_integrated_boundary_polynomial_witness_probe"]
        cls.witness = cls.probe["track_B_triaxial_ellipsoid_integrated_boundary_polynomial_witness"]
        cls.ellipsoid = cls.witness["ellipsoid_witness"]

    def test_packet_identity_and_outputs(self) -> None:
        self.assertEqual(self.payload["schema_version"], "p2349_s1299_v1")
        self.assertEqual(self.payload["packet_id"], "P2349")
        self.assertEqual(self.payload["stage_id"], "S1299")
        self.assertTrue(OUT.exists())
        self.assertTrue(MD.exists())

    def test_axis_pole_samples_and_integrated_values(self) -> None:
        self.assertEqual(self.witness["ledger_coefficient_b_GB"], "13152087*log(2)/(320000000*pi**2)")
        self.assertEqual(self.witness["per_euler_pairing"], "13152087*log(2)/10000000")
        self.assertEqual(self.ellipsoid["axes"], ["1", "2", "3", "5"])
        self.assertEqual(self.ellipsoid["distinct_sigma3_count"], 4)
        self.assertEqual(self.ellipsoid["distinct_boundary_density_count"], 4)
        rows = self.ellipsoid["axis_pole_rows"]
        self.assertEqual(len(rows), 4)
        self.assertEqual(rows[0]["principal_curvatures"], ["1/4", "1/9", "1/25"])
        self.assertEqual(rows[0]["sigma3"], "1/900")
        self.assertEqual(rows[3]["principal_curvatures"], ["5", "5/4", "5/9"])
        self.assertEqual(rows[3]["sigma3"], "125/36")
        self.assertEqual(self.ellipsoid["nonconstant_density_witness_x4_minus_x1"], "12496/225")
        self.assertEqual(self.ellipsoid["gauss_map_degree"], "1")
        self.assertEqual(self.ellipsoid["integral_sigma3_via_gauss_map"], "2*pi**2")
        self.assertEqual(self.ellipsoid["integrated_boundary_functional_16_integral_sigma3"], "32*pi**2")
        self.assertEqual(self.ellipsoid["integrated_boundary_residual"], "0")
        self.assertEqual(self.ellipsoid["normalized_track_B_pairing"], "13152087*log(2)/10000000")
        self.assertEqual(self.ellipsoid["normalized_pairing_residual"], "0")

    def test_dependencies_gatekeepers_and_fingerprint(self) -> None:
        deps = self.probe["current_export_dependencies"]
        self.assertTrue(deps["p2345_boundary_value_replayed"])
        self.assertTrue(deps["p2348_flat_integral_replay_replayed"])
        for key, value in self.payload["gatekeeper_checks"].items():
            self.assertTrue(value, key)
        theorem = self.probe["theorem_export"]
        self.assertEqual(self.probe["theorem_fingerprint_sha256"], sha256_json(theorem))
        self.assertIn("integrated arbitrary-boundary theorem", theorem["not_licensed"])
        self.assertIn("selector premise or QW-2191 selector discharge", theorem["not_licensed"])
        self.assertIn("ToE closure", theorem["not_licensed"])


if __name__ == "__main__":
    unittest.main()
