from __future__ import annotations

import hashlib
import json
import subprocess
import sys
import unittest
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parent
SCRIPT = ROOT / "p2360_s1310_strict_track_b_smooth_interface_chern_polynomial_sign_reversal_lemma.py"
OUT = ROOT / "generated" / "p2360_s1310_strict_track_b_smooth_interface_chern_polynomial_sign_reversal_lemma.json"
MD = ROOT / "generated" / "p2360_s1310_strict_track_b_smooth_interface_chern_polynomial_sign_reversal_lemma.md"


def sha256_json(obj: Any) -> str:
    return hashlib.sha256(
        json.dumps(obj, sort_keys=True, ensure_ascii=False, separators=(",", ":")).encode("utf-8")
    ).hexdigest()


class P2360SmoothInterfaceChernPolynomialSignReversalTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        subprocess.run([sys.executable, str(SCRIPT)], check=True, cwd=ROOT.parent)
        cls.payload = json.loads(OUT.read_text(encoding="utf-8"))
        cls.probe = cls.payload["strict_track_b_smooth_interface_chern_polynomial_sign_reversal_lemma_probe"]
        cls.lemma = cls.probe["track_B_smooth_interface_chern_polynomial_sign_reversal_lemma"]

    def test_packet_identity_and_outputs(self) -> None:
        self.assertEqual(self.payload["schema_version"], "p2360_s1310_v1")
        self.assertEqual(self.payload["packet_id"], "P2360")
        self.assertEqual(self.payload["stage_id"], "S1310")
        self.assertTrue(OUT.exists())
        self.assertTrue(MD.exists())

    def test_symbolic_sign_reversal(self) -> None:
        self.assertEqual(self.lemma["chern_boundary_polynomial"], "8*(K*k1 + K*k2 + K*k3 + 2*k1*k2*k3)")
        self.assertEqual(self.lemma["reversed_normal_polynomial"], "-8*(K*k1 + K*k2 + K*k3 + 2*k1*k2*k3)")
        self.assertEqual(self.lemma["sign_reversal_residual"], "0")
        self.assertEqual(self.lemma["interface_pair_sum"], "0")
        self.assertEqual(self.lemma["integrated_interface_residual_symbolic"], "0")
        self.assertEqual(self.lemma["pairing_interface_residual_symbolic"], "0")
        self.assertEqual(self.lemma["sigma1_reversal_residual"], "0")
        self.assertEqual(self.lemma["sigma2_even_reversal_residual_recorded_not_used"], "0")
        self.assertEqual(self.lemma["sigma3_reversal_residual"], "0")
        self.assertEqual(self.lemma["p2348_polynomial_residual"], "0")

    def test_samples_and_upstream_replays(self) -> None:
        rows = {row["sample_id"]: row for row in self.lemma["sample_rows"]}
        self.assertEqual(len(rows), 3)
        self.assertTrue(self.lemma["all_sample_sums_zero"])
        self.assertEqual(rows["p2348_nonsymmetric_sample_K5_123"]["value"], "336")
        self.assertEqual(rows["p2348_nonsymmetric_sample_K5_123"]["reversed_normal_value"], "-336")
        self.assertEqual(rows["mixed_curvature_sample_K_minus2_347"]["value"], "1120")
        self.assertEqual(rows["zero_sigma3_nonzero_sigma1_sample"]["value"], "616")
        for row in rows.values():
            self.assertEqual(row["interface_sum"], "0", row["sample_id"])
        self.assertEqual(self.lemma["p2358_symbolic_interface_residual_replayed"], "0")
        self.assertTrue(self.lemma["p2358_all_case_residuals_replayed"])
        self.assertEqual(self.lemma["p2359_interface_cancellation_replayed"], "0")
        self.assertEqual(self.lemma["p2359_gluing_consistency_residual_replayed"], "0")
        self.assertIn("O5_regularization_corners_and_gluing", self.lemma["p2353_minimal_cut_replayed"])
        self.assertTrue(self.lemma["o5_smooth_interface_polynomial_partially_attacked"])

    def test_dependencies_gatekeepers_and_fingerprint(self) -> None:
        for key, value in self.probe["current_export_dependencies"].items():
            self.assertTrue(value, key)
        for key, value in self.payload["gatekeeper_checks"].items():
            self.assertTrue(value, key)
        theorem = self.probe["theorem_export"]
        self.assertEqual(self.probe["theorem_fingerprint_sha256"], sha256_json(theorem))
        self.assertIn("arbitrary-boundary theorem", theorem["not_licensed"])
        self.assertIn("corner contribution theorem", theorem["not_licensed"])
        self.assertIn("selector premise or QW-2191 selector discharge", theorem["not_licensed"])
        self.assertIn("ToE closure", theorem["not_licensed"])


if __name__ == "__main__":
    unittest.main()
