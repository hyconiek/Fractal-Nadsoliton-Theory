from __future__ import annotations

import hashlib
import json
import subprocess
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent
G = ROOT / "generated"


def sha256_json(payload: object) -> str:
    blob = json.dumps(payload, sort_keys=True, separators=(",", ":"))
    return hashlib.sha256(blob.encode("utf-8")).hexdigest()


class TestP2312S1262StrictLapseShearSourceNormalizationAttemptProbe(unittest.TestCase):
    def test_packet(self) -> None:
        subprocess.run(
            [sys.executable, str(ROOT / "p2312_s1262_strict_lapse_shear_source_normalization_attempt_probe.py")],
            check=True,
        )
        data = json.loads(
            (G / "p2312_s1262_strict_lapse_shear_source_normalization_attempt_probe.json").read_text(encoding="utf-8")
        )
        self.assertEqual(data["schema_version"], "p2312_s1262_v1")
        self.assertEqual(data["packet_id"], "P2312")
        probe = data["strict_lapse_shear_source_normalization_attempt_probe"]
        candidate = probe["candidate_source"]
        self.assertEqual(candidate["source_id"], "P1985_STRICT_NON_GB_LAPSE_RESIDUAL")
        self.assertFalse(candidate["strict_admissible"])
        witnesses = probe["orientation_witnesses"]
        signs = {row["source_sign"] for row in witnesses}
        self.assertIn("positive", signs)
        self.assertIn("negative", signs)
        self.assertEqual({row["shear_energy_Q"] for row in witnesses}, {1.0})
        attempts = probe["normalization_attempts"]
        self.assertFalse(any(row["strict_admissible"] for row in attempts))
        verdict = probe["source_normalization_verdict"]
        self.assertFalse(verdict["new_lapse_shear_source_found"])
        self.assertTrue(verdict["direct_non_gb_lapse_source_sign_indefinite"])
        self.assertTrue(verdict["same_shear_energy_allows_opposite_source_signs"])
        self.assertTrue(verdict["self_energy_numeric_condition_passes"])
        self.assertFalse(verdict["admissible_for_provider_lift_export"])
        self.assertFalse(verdict["admissible_for_g1_g3_update"])
        task3 = probe["task3_g1_g3_update"]
        self.assertEqual(task3["G1_reduction_certainty"], "OPEN_UNCHANGED")
        self.assertEqual(task3["G2_nonlinear_trajectory_realism"], "CLOSED_FROM_P2301_UNCHANGED")
        self.assertEqual(task3["G3_operational_policy_rule"], "OPEN_UNCHANGED")
        theorem = probe["theorem_export"]
        self.assertTrue(theorem["proof_bits"]["sign_indefinite"])
        self.assertTrue(theorem["proof_bits"]["same_shear_energy_sign_flip"])
        self.assertEqual(sha256_json(theorem), probe["theorem_fingerprint_sha256"])
        g = data["gatekeeper_checks"]
        self.assertTrue(g["p1985_non_gb_lapse_operator_loaded"])
        self.assertTrue(g["p1986_non_gb_decomposition_loaded"])
        self.assertTrue(g["p2300_coefficients_loaded"])
        self.assertTrue(g["p2310_self_energy_numeric_condition_loaded"])
        self.assertTrue(g["p2311_bridge_blocker_loaded"])
        self.assertTrue(g["explicit_positive_and_negative_lapse_source_witnesses"])
        self.assertTrue(g["same_shear_energy_sign_flip_witnessed"])
        self.assertTrue(g["all_normalization_attempts_quarantined"])
        self.assertTrue(g["no_new_lapse_shear_source_exported"])
        self.assertTrue(g["strict_g1_g3_not_updated"])
        self.assertTrue(g["no_selector_premise_added"])
        self.assertTrue(g["no_qw2191_discharge_claimed"])
        self.assertTrue(g["no_legacy_kernel_role_transfer"])
        self.assertTrue(g["no_toe_closure_claimed"])


if __name__ == "__main__":
    unittest.main()
