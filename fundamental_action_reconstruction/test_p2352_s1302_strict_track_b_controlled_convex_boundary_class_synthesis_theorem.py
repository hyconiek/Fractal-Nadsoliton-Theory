from __future__ import annotations

import hashlib
import json
import subprocess
import sys
import unittest
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parent
SCRIPT = ROOT / "p2352_s1302_strict_track_b_controlled_convex_boundary_class_synthesis_theorem.py"
OUT = ROOT / "generated" / "p2352_s1302_strict_track_b_controlled_convex_boundary_class_synthesis_theorem.json"
MD = ROOT / "generated" / "p2352_s1302_strict_track_b_controlled_convex_boundary_class_synthesis_theorem.md"


def sha256_json(obj: Any) -> str:
    return hashlib.sha256(
        json.dumps(obj, sort_keys=True, ensure_ascii=False, separators=(",", ":")).encode("utf-8")
    ).hexdigest()


class P2352ControlledConvexBoundaryClassSynthesisTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        subprocess.run([sys.executable, str(SCRIPT)], check=True, cwd=ROOT.parent)
        cls.payload = json.loads(OUT.read_text(encoding="utf-8"))
        cls.probe = cls.payload["strict_track_b_controlled_convex_boundary_class_synthesis_theorem_probe"]
        cls.synthesis = cls.probe["track_B_controlled_convex_boundary_class_synthesis"]

    def test_packet_identity_and_outputs(self) -> None:
        self.assertEqual(self.payload["schema_version"], "p2352_s1302_v1")
        self.assertEqual(self.payload["packet_id"], "P2352")
        self.assertEqual(self.payload["stage_id"], "S1302")
        self.assertTrue(OUT.exists())
        self.assertTrue(MD.exists())

    def test_controlled_rows_and_residual_vectors(self) -> None:
        self.assertEqual(self.synthesis["ledger_coefficient_b_GB"], "13152087*log(2)/(320000000*pi**2)")
        self.assertEqual(self.synthesis["target_boundary_functional"], "32*pi**2")
        self.assertEqual(self.synthesis["target_pairing_per_euler_number"], "13152087*log(2)/10000000")
        rows = self.synthesis["class_rows"]
        self.assertEqual(len(rows), 5)
        self.assertEqual([row["packet"] for row in rows], ["P2343", "P2344", "P2349", "P2350", "P2351"])
        self.assertTrue(all(row["boundary_functional"] == "32*pi**2" for row in rows))
        self.assertEqual(self.synthesis["boundary_residual_vector"], ["0", "0", "0", "0", "0"])
        self.assertEqual(self.synthesis["pairing_residual_vector"], ["0", "0", "0", "0", "0"])
        self.assertTrue(self.synthesis["all_boundary_residuals_zero"])
        self.assertTrue(self.synthesis["all_pairing_residuals_zero"])

    def test_feature_matrix_dependencies_gatekeepers_and_fingerprint(self) -> None:
        self.assertEqual(self.synthesis["feature_rank"], 4)
        self.assertEqual(self.synthesis["nonround_or_nonconstant_count"], 4)
        self.assertEqual(self.synthesis["global_convexity_certificate_count"], 4)
        self.assertEqual(self.synthesis["local_chern_bridge_count"], 3)
        self.assertEqual(self.synthesis["coverage_score_raw"], "11")
        self.assertEqual(self.synthesis["coverage_score_normalized"], "11/15")
        self.assertEqual(self.synthesis["p2345_convex_class_boundary_residual"], "0")
        self.assertEqual(self.synthesis["p2348_chern_polynomial_replay_residual"], "0")
        for key, value in self.probe["current_export_dependencies"].items():
            self.assertTrue(value, key)
        for key, value in self.payload["gatekeeper_checks"].items():
            self.assertTrue(value, key)
        theorem = self.probe["theorem_export"]
        self.assertEqual(self.probe["theorem_fingerprint_sha256"], sha256_json(theorem))
        self.assertIn("arbitrary-boundary theorem", theorem["not_licensed"])
        self.assertIn("selector premise or QW-2191 selector discharge", theorem["not_licensed"])
        self.assertIn("ToE closure", theorem["not_licensed"])


if __name__ == "__main__":
    unittest.main()
