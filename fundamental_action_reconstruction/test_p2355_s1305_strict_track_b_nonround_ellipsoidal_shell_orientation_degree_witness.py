from __future__ import annotations

import hashlib
import json
import subprocess
import sys
import unittest
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parent
SCRIPT = ROOT / "p2355_s1305_strict_track_b_nonround_ellipsoidal_shell_orientation_degree_witness.py"
OUT = ROOT / "generated" / "p2355_s1305_strict_track_b_nonround_ellipsoidal_shell_orientation_degree_witness.json"
MD = ROOT / "generated" / "p2355_s1305_strict_track_b_nonround_ellipsoidal_shell_orientation_degree_witness.md"


def sha256_json(obj: Any) -> str:
    return hashlib.sha256(
        json.dumps(obj, sort_keys=True, ensure_ascii=False, separators=(",", ":")).encode("utf-8")
    ).hexdigest()


class P2355NonroundEllipsoidalShellOrientationDegreeTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        subprocess.run([sys.executable, str(SCRIPT)], check=True, cwd=ROOT.parent)
        cls.payload = json.loads(OUT.read_text(encoding="utf-8"))
        cls.probe = cls.payload["strict_track_b_nonround_ellipsoidal_shell_orientation_degree_witness_probe"]
        cls.witness = cls.probe["track_B_nonround_ellipsoidal_shell_orientation_degree_witness"]

    def test_packet_identity_and_outputs(self) -> None:
        self.assertEqual(self.payload["schema_version"], "p2355_s1305_v1")
        self.assertEqual(self.payload["packet_id"], "P2355")
        self.assertEqual(self.payload["stage_id"], "S1305")
        self.assertTrue(OUT.exists())
        self.assertTrue(MD.exists())

    def test_axis_samples_are_nonround_and_signed(self) -> None:
        self.assertEqual(self.witness["strict_coordinate_axis_containment_checks"], [True, True, True, True])
        self.assertTrue(self.witness["both_components_have_nonconstant_abs_sigma3"])
        outer, inner = self.witness["components"]
        self.assertEqual(outer["axes"], ["3", "4", "6", "9"])
        self.assertEqual(inner["axes"], ["1", "2", "3", "5"])
        self.assertEqual(outer["signed_sigma3_values"], ["1/1728", "16/6561", "1/54", "9/64"])
        self.assertEqual(inner["signed_sigma3_values"], ["-1/900", "-8/225", "-27/100", "-125/36"])
        self.assertTrue(outer["nonconstant_abs_sigma3"])
        self.assertTrue(inner["nonconstant_abs_sigma3"])
        self.assertEqual(outer["axis_pole_samples"][0]["boundary_density_16_sigma3"], "1/108")
        self.assertEqual(inner["axis_pole_samples"][3]["boundary_density_16_sigma3"], "-500/9")

    def test_shell_cancellation_residuals_and_replays(self) -> None:
        self.assertEqual(self.witness["target_boundary_per_degree"], "32*pi**2")
        self.assertEqual(self.witness["target_pairing_per_degree"], "13152087*log(2)/10000000")
        self.assertEqual(self.witness["component_gauss_degree_sum"], 0)
        self.assertEqual(self.witness["component_boundary_functional_sum"], "0")
        self.assertEqual(self.witness["component_pairing_sum"], "0")
        self.assertEqual(self.witness["shell_boundary_residual"], "0")
        self.assertEqual(self.witness["shell_pairing_residual"], "0")
        self.assertEqual(self.witness["degree_boundary_residual"], "0")
        self.assertEqual(self.witness["degree_pairing_residual"], "0")
        self.assertEqual(self.witness["p2345_convex_boundary_residual"], "0")
        self.assertEqual(self.witness["p2348_chern_replay_residual"], "0")
        self.assertEqual(self.witness["p2349_triaxial_replay_residual"], "0")
        self.assertEqual(self.witness["p2354_shell_boundary_residual_replayed"], "0")
        self.assertIn("O4_nonconvex_degree_and_orientation_accounting", self.witness["p2353_minimal_cut_replayed"])

    def test_dependencies_gatekeepers_and_fingerprint(self) -> None:
        for key, value in self.probe["current_export_dependencies"].items():
            self.assertTrue(value, key)
        for key, value in self.payload["gatekeeper_checks"].items():
            self.assertTrue(value, key)
        theorem = self.probe["theorem_export"]
        self.assertEqual(self.probe["theorem_fingerprint_sha256"], sha256_json(theorem))
        self.assertIn("arbitrary-boundary theorem", theorem["not_licensed"])
        self.assertIn("general nonconvex boundary theorem", theorem["not_licensed"])
        self.assertIn("selector premise or QW-2191 selector discharge", theorem["not_licensed"])
        self.assertIn("ToE closure", theorem["not_licensed"])


if __name__ == "__main__":
    unittest.main()
