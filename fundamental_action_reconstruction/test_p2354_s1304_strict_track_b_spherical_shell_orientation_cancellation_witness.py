from __future__ import annotations

import hashlib
import json
import subprocess
import sys
import unittest
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parent
SCRIPT = ROOT / "p2354_s1304_strict_track_b_spherical_shell_orientation_cancellation_witness.py"
OUT = ROOT / "generated" / "p2354_s1304_strict_track_b_spherical_shell_orientation_cancellation_witness.json"
MD = ROOT / "generated" / "p2354_s1304_strict_track_b_spherical_shell_orientation_cancellation_witness.md"


def sha256_json(obj: Any) -> str:
    return hashlib.sha256(
        json.dumps(obj, sort_keys=True, ensure_ascii=False, separators=(",", ":")).encode("utf-8")
    ).hexdigest()


class P2354SphericalShellOrientationCancellationTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        subprocess.run([sys.executable, str(SCRIPT)], check=True, cwd=ROOT.parent)
        cls.payload = json.loads(OUT.read_text(encoding="utf-8"))
        cls.probe = cls.payload["strict_track_b_spherical_shell_orientation_cancellation_witness_probe"]
        cls.witness = cls.probe["track_B_spherical_shell_orientation_cancellation_witness"]

    def test_packet_identity_and_outputs(self) -> None:
        self.assertEqual(self.payload["schema_version"], "p2354_s1304_v1")
        self.assertEqual(self.payload["packet_id"], "P2354")
        self.assertEqual(self.payload["stage_id"], "S1304")
        self.assertTrue(OUT.exists())
        self.assertTrue(MD.exists())

    def test_component_orientation_values(self) -> None:
        components = self.witness["components"]
        self.assertEqual(len(components), 2)
        outer, inner = components
        self.assertEqual(outer["principal_curvatures"], ["1/3", "1/3", "1/3"])
        self.assertEqual(outer["sigma3"], "1/27")
        self.assertEqual(outer["area"], "54*pi**2")
        self.assertEqual(outer["integral_sigma3_dA"], "2*pi**2")
        self.assertEqual(outer["boundary_functional_16_integral_sigma3"], "32*pi**2")
        self.assertEqual(outer["gauss_map_degree"], 1)
        self.assertEqual(inner["principal_curvatures"], ["-2", "-2", "-2"])
        self.assertEqual(inner["sigma3"], "-8")
        self.assertEqual(inner["area"], "pi**2/4")
        self.assertEqual(inner["integral_sigma3_dA"], "-2*pi**2")
        self.assertEqual(inner["boundary_functional_16_integral_sigma3"], "-32*pi**2")
        self.assertEqual(inner["gauss_map_degree"], -1)

    def test_shell_cancellation_residuals_and_o4_replay(self) -> None:
        self.assertEqual(self.witness["target_boundary_per_degree"], "32*pi**2")
        self.assertEqual(self.witness["target_pairing_per_degree"], "13152087*log(2)/10000000")
        self.assertEqual(self.witness["component_gauss_degree_sum"], 0)
        self.assertEqual(self.witness["component_boundary_functional_sum"], "0")
        self.assertEqual(self.witness["component_pairing_sum"], "0")
        self.assertEqual(self.witness["shell_boundary_residual"], "0")
        self.assertEqual(self.witness["shell_pairing_residual"], "0")
        self.assertEqual(self.witness["degree_boundary_residual"], "0")
        self.assertEqual(self.witness["degree_pairing_residual"], "0")
        self.assertTrue(self.witness["o4_cut_partially_attacked"])
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
