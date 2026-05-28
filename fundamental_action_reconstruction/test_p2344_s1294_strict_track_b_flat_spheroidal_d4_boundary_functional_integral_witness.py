from __future__ import annotations

import hashlib
import json
import subprocess
import sys
import unittest
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parent
SCRIPT = ROOT / "p2344_s1294_strict_track_b_flat_spheroidal_d4_boundary_functional_integral_witness.py"
OUT = ROOT / "generated" / "p2344_s1294_strict_track_b_flat_spheroidal_d4_boundary_functional_integral_witness.json"
MD = ROOT / "generated" / "p2344_s1294_strict_track_b_flat_spheroidal_d4_boundary_functional_integral_witness.md"


def sha256_json(obj: Any) -> str:
    return hashlib.sha256(
        json.dumps(obj, sort_keys=True, ensure_ascii=False, separators=(",", ":")).encode("utf-8")
    ).hexdigest()


class P2344TrackBFlatSpheroidalD4BoundaryFunctionalTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        subprocess.run([sys.executable, str(SCRIPT)], check=True, cwd=ROOT.parent)
        cls.payload = json.loads(OUT.read_text(encoding="utf-8"))
        cls.probe = cls.payload["strict_track_b_flat_spheroidal_d4_boundary_functional_integral_witness_probe"]
        cls.witness = cls.probe["track_B_flat_spheroidal_D4_boundary_functional_integral_witness"]
        cls.spheroidal = cls.witness["spheroidal_boundary_class"]

    def test_packet_identity_and_outputs(self) -> None:
        self.assertEqual(self.payload["schema_version"], "p2344_s1294_v1")
        self.assertEqual(self.payload["packet_id"], "P2344")
        self.assertEqual(self.payload["stage_id"], "S1294")
        self.assertTrue(OUT.exists())
        self.assertTrue(MD.exists())

    def test_spheroidal_boundary_integral_values(self) -> None:
        self.assertEqual(self.witness["ledger_coefficient_b_GB"], "13152087*log(2)/(320000000*pi**2)")
        self.assertEqual(self.witness["per_euler_pairing"], "13152087*log(2)/10000000")
        self.assertEqual(self.spheroidal["nonround_subclass_condition"], "q != 1")
        self.assertEqual(self.spheroidal["bulk_gb_fixture_value"], "0")
        self.assertEqual(
            self.spheroidal["gauss_map_jacobian_density"],
            "q**3/(q**2*sin(theta)**2 + cos(theta)**2)**2",
        )
        self.assertEqual(self.spheroidal["half_interval_t_integrand"], "q**4*t**2/(q**2*t**2 + 1)**2")
        self.assertEqual(self.spheroidal["half_interval_integral_value"], "pi*q/4")
        self.assertEqual(self.spheroidal["sigma3_integral_over_boundary"], "2*pi**2")
        self.assertEqual(self.spheroidal["sigma3_integral_residual"], "0")
        self.assertEqual(self.spheroidal["boundary_functional_B_spheroidal_q"], "32*pi**2")
        self.assertEqual(self.spheroidal["boundary_functional_residual"], "0")
        self.assertEqual(self.spheroidal["q_derivative_residual"], "0")
        self.assertEqual(self.spheroidal["round_limit_value_q1"], "32*pi**2")
        self.assertEqual(self.spheroidal["sample_q"], "2")
        self.assertEqual(self.spheroidal["sample_q2_boundary_functional"], "32*pi**2")
        self.assertEqual(self.spheroidal["sample_q2_density_pole"], "8")
        self.assertEqual(self.spheroidal["sample_q2_density_equator"], "1/2")
        self.assertEqual(self.spheroidal["sample_q2_nonconstant_density_difference"], "15/2")
        self.assertEqual(self.spheroidal["normalized_track_B_pairing"], "13152087*log(2)/10000000")
        self.assertEqual(self.spheroidal["normalized_pairing_residual"], "0")

    def test_dependencies_gatekeepers_and_fingerprint(self) -> None:
        deps = self.probe["current_export_dependencies"]
        self.assertTrue(deps["p2341_boundary_value_rederived"])
        self.assertTrue(deps["p2341_pairing_rederived"])
        self.assertTrue(deps["p2343_round_limit_rederived"])
        for key, value in self.payload["gatekeeper_checks"].items():
            self.assertTrue(value, key)
        theorem = self.probe["theorem_export"]
        self.assertEqual(self.probe["theorem_fingerprint_sha256"], sha256_json(theorem))
        self.assertIn("boundary functional for arbitrary non-spheroidal or non-convex boundaries", theorem["not_licensed"])
        self.assertIn("selector premise or QW-2191 selector discharge", theorem["not_licensed"])
        self.assertIn("ToE closure", theorem["not_licensed"])


if __name__ == "__main__":
    unittest.main()
