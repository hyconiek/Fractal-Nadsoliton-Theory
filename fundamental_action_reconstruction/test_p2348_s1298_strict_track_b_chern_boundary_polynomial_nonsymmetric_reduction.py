from __future__ import annotations

import hashlib
import json
import subprocess
import sys
import unittest
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parent
SCRIPT = ROOT / "p2348_s1298_strict_track_b_chern_boundary_polynomial_nonsymmetric_reduction.py"
OUT = ROOT / "generated" / "p2348_s1298_strict_track_b_chern_boundary_polynomial_nonsymmetric_reduction.json"
MD = ROOT / "generated" / "p2348_s1298_strict_track_b_chern_boundary_polynomial_nonsymmetric_reduction.md"


def sha256_json(obj: Any) -> str:
    return hashlib.sha256(
        json.dumps(obj, sort_keys=True, ensure_ascii=False, separators=(",", ":")).encode("utf-8")
    ).hexdigest()


class P2348ChernBoundaryPolynomialNonsymmetricTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        subprocess.run([sys.executable, str(SCRIPT)], check=True, cwd=ROOT.parent)
        cls.payload = json.loads(OUT.read_text(encoding="utf-8"))
        cls.probe = cls.payload["strict_track_b_chern_boundary_polynomial_nonsymmetric_reduction_probe"]
        cls.witness = cls.probe["track_B_chern_boundary_polynomial_nonsymmetric_reduction"]
        cls.poly = cls.witness["boundary_polynomial"]

    def test_packet_identity_and_outputs(self) -> None:
        self.assertEqual(self.payload["schema_version"], "p2348_s1298_v1")
        self.assertEqual(self.payload["packet_id"], "P2348")
        self.assertEqual(self.payload["stage_id"], "S1298")
        self.assertTrue(OUT.exists())
        self.assertTrue(MD.exists())

    def test_polynomial_reductions_and_sample(self) -> None:
        self.assertEqual(self.witness["ledger_coefficient_b_GB"], "13152087*log(2)/(320000000*pi**2)")
        self.assertEqual(self.witness["per_euler_pairing"], "13152087*log(2)/10000000")
        self.assertEqual(self.poly["sigma1"], "k1 + k2 + k3")
        self.assertEqual(self.poly["sigma2_recorded_not_used"], "k1*k2 + k1*k3 + k2*k3")
        self.assertEqual(self.poly["sigma3"], "k1*k2*k3")
        self.assertEqual(self.poly["boundary_polynomial"], "8*(K*k1 + K*k2 + K*k3 + 2*k1*k2*k3)")
        self.assertEqual(self.poly["flat_K0_reduction"], "16*k1*k2*k3")
        self.assertEqual(self.poly["flat_K0_residual"], "0")
        self.assertEqual(self.poly["spherical_cap_equal_curvature_reduction"], "8*k*(2*k**2 + 3)")
        self.assertEqual(self.poly["spherical_cap_residual"], "0")
        self.assertTrue(self.poly["all_permutation_residuals_zero"])
        sample = self.poly["nonsymmetric_sample"]
        self.assertEqual(sample["sigma1"], "6")
        self.assertEqual(sample["sigma2_recorded_not_used"], "11")
        self.assertEqual(sample["sigma3"], "6")
        self.assertEqual(sample["boundary_polynomial_value"], "336")
        self.assertEqual(sample["residual"], "0")

    def test_flat_integral_replay_dependencies_gatekeepers_and_fingerprint(self) -> None:
        replay = self.poly["flat_convex_integral_replay"]
        self.assertEqual(replay["substitute_integral_sigma3_equals_2pi2"], "32*pi**2")
        self.assertEqual(replay["flat_convex_residual"], "0")
        self.assertEqual(replay["flat_convex_pairing"], "13152087*log(2)/10000000")
        self.assertEqual(replay["flat_convex_pairing_residual"], "0")
        deps = self.probe["current_export_dependencies"]
        self.assertTrue(deps["p2343_flat_value_replayed"])
        self.assertTrue(deps["p2345_convex_value_replayed"])
        for key, value in self.payload["gatekeeper_checks"].items():
            self.assertTrue(value, key)
        theorem = self.probe["theorem_export"]
        self.assertEqual(self.probe["theorem_fingerprint_sha256"], sha256_json(theorem))
        self.assertIn("integrated arbitrary-boundary theorem", theorem["not_licensed"])
        self.assertIn("selector premise or QW-2191 selector discharge", theorem["not_licensed"])
        self.assertIn("ToE closure", theorem["not_licensed"])


if __name__ == "__main__":
    unittest.main()
