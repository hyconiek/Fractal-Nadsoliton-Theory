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


class TestP2302S1252StrictTask3ProviderLiftPolicyLockCandidateProbe(unittest.TestCase):
    def test_packet(self) -> None:
        subprocess.run(
            [sys.executable, str(ROOT / "p2302_s1252_strict_task3_provider_lift_policy_lock_candidate_probe.py")],
            check=True,
        )
        data = json.loads(
            (G / "p2302_s1252_strict_task3_provider_lift_policy_lock_candidate_probe.json").read_text(encoding="utf-8")
        )
        self.assertEqual(data["schema_version"], "p2302_s1252_v1")
        probe = data["strict_task3_provider_lift_policy_lock_candidate_probe"]
        candidate = probe["policy_lock_candidate"]
        self.assertEqual(candidate["status"], "CONDITIONAL_CANDIDATE_REQUIRES_PROVIDER_TO_MARGIN_BRIDGE")
        self.assertGreater(candidate["provider_lift_per_step"], 0.0)
        self.assertLess(candidate["provider_lift_per_step"], 0.01)
        self.assertTrue(probe["conditional_fresh_replay"]["summary"]["all_rows_meet_target"])
        self.assertEqual(probe["strict_closure_status"], "HELD_OPEN_UNTIL_PROVIDER_TO_MARGIN_BRIDGE_PROVEN")
        self.assertEqual(sha256_json(probe["theorem_export"]), probe["theorem_fingerprint_sha256"])
        g = data["gatekeeper_checks"]
        self.assertTrue(g["alpha_geo_strict_source_loaded"])
        self.assertTrue(g["alpha_geo_is_four_ln2_not_legacy_import"])
        self.assertTrue(g["p2300_canonical_coefficients_loaded"])
        self.assertTrue(g["p2301_only_g2_closed_input"])
        self.assertTrue(g["provider_lift_candidate_found"])
        self.assertTrue(g["conditional_replay_all_rows_meet_target"])
        self.assertTrue(g["baseline_p2280_had_no_feasible_lock"])
        self.assertTrue(g["provider_to_margin_bridge_still_open"])
        self.assertTrue(g["strict_task3_closure_not_claimed"])
        self.assertTrue(g["no_qw2191_discharge_claimed"])
        self.assertTrue(g["no_selector_closure_claimed"])
        self.assertTrue(g["no_legacy_kernel_role_transfer"])
        self.assertTrue(g["no_toe_closure_claimed"])


if __name__ == "__main__":
    unittest.main()
