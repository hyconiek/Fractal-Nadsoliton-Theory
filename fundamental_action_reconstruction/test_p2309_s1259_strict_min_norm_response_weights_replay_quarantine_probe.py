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


class TestP2309S1259StrictMinNormResponseWeightsReplayQuarantineProbe(unittest.TestCase):
    def test_packet(self) -> None:
        subprocess.run(
            [sys.executable, str(ROOT / "p2309_s1259_strict_min_norm_response_weights_replay_quarantine_probe.py")],
            check=True,
        )
        data = json.loads(
            (G / "p2309_s1259_strict_min_norm_response_weights_replay_quarantine_probe.json").read_text(encoding="utf-8")
        )
        self.assertEqual(data["schema_version"], "p2309_s1259_v1")
        probe = data["strict_min_norm_response_weights_replay_quarantine_probe"]
        attempt = probe["variational_identity_attempt"]
        self.assertTrue(attempt["weights_exported"])
        self.assertFalse(attempt["strict_admissible"])
        self.assertAlmostEqual(attempt["induced_lift"], attempt["required_lift"], places=12)
        self.assertEqual(len(attempt["weights"]), 10)
        replay = probe["target_calibrated_replay"]["summary"]
        self.assertTrue(replay["all_rows_meet_target"])
        self.assertGreater(replay["worst_margin_to_target"], 0.0)
        verdict = probe["strict_admissibility_verdict"]
        self.assertTrue(verdict["lambda_equals_R_of_c_exported"])
        self.assertFalse(verdict["lambda_equals_R_of_c_strictly_derived"])
        self.assertTrue(verdict["replay_passes"])
        self.assertFalse(verdict["admissible_for_g1_g3_update"])
        self.assertEqual(probe["task3_g1_g3_update"]["G1_reduction_certainty"], "OPEN_UNCHANGED")
        self.assertEqual(probe["task3_g1_g3_update"]["G3_operational_policy_rule"], "OPEN_UNCHANGED")
        self.assertEqual(sha256_json(probe["theorem_export"]), probe["theorem_fingerprint_sha256"])
        g = data["gatekeeper_checks"]
        self.assertTrue(g["alpha_geo_strict_source_loaded"])
        self.assertTrue(g["alpha_geo_is_four_ln2_not_legacy_import"])
        self.assertTrue(g["p2300_coefficients_loaded"])
        self.assertTrue(g["p2302_required_lift_loaded"])
        self.assertTrue(g["p2308_current_class_nonexistence_loaded"])
        self.assertTrue(g["min_norm_weights_exported"])
        self.assertTrue(g["induced_lift_matches_required_lift"])
        self.assertTrue(g["replay_passes_if_target_calibrated_weights_are_admitted"])
        self.assertTrue(g["weights_are_target_calibrated_not_strictly_derived"])
        self.assertTrue(g["strict_admissible_response_bridge_not_claimed"])
        self.assertTrue(g["strict_g1_g3_not_updated"])
        self.assertTrue(g["no_qw2191_discharge_claimed"])
        self.assertTrue(g["no_selector_closure_claimed"])
        self.assertTrue(g["no_legacy_kernel_role_transfer"])
        self.assertTrue(g["no_toe_closure_claimed"])


if __name__ == "__main__":
    unittest.main()
