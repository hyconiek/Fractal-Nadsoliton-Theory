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


class TestP2304S1254StrictDirectFloorPolicyLockReplayProbe(unittest.TestCase):
    def test_packet(self) -> None:
        subprocess.run(
            [sys.executable, str(ROOT / "p2304_s1254_strict_direct_floor_policy_lock_replay_probe.py")],
            check=True,
        )
        data = json.loads(
            (G / "p2304_s1254_strict_direct_floor_policy_lock_replay_probe.json").read_text(encoding="utf-8")
        )
        self.assertEqual(data["schema_version"], "p2304_s1254_v1")
        probe = data["strict_direct_floor_policy_lock_replay_probe"]
        comparison = probe["lift_comparison"]
        self.assertLess(comparison["p2303_strict_direct_floor_lift"], comparison["p2302_required_lift"])
        summary = probe["direct_floor_fresh_replay"]["summary"]
        self.assertFalse(summary["all_rows_meet_target"])
        self.assertLess(summary["worst_margin_to_target"], 0.0)
        self.assertEqual([row["status"] for row in probe["direct_floor_gap_rows"]], ["OPEN", "CLOSED_FROM_P2301", "OPEN"])
        self.assertEqual(probe["strict_closure_status"], "HELD_OPEN_DIRECT_FLOOR_ROUTE_REFUTED")
        self.assertEqual(sha256_json(probe["theorem_export"]), probe["theorem_fingerprint_sha256"])
        g = data["gatekeeper_checks"]
        self.assertTrue(g["alpha_geo_strict_source_loaded"])
        self.assertTrue(g["alpha_geo_is_four_ln2_not_legacy_import"])
        self.assertTrue(g["p2303_direct_floor_loaded"])
        self.assertTrue(g["p2302_required_lift_loaded"])
        self.assertTrue(g["direct_floor_below_required_lift"])
        self.assertTrue(g["direct_floor_replay_rows_exported"])
        self.assertTrue(g["direct_floor_replay_fails_g1_g3"])
        self.assertTrue(g["direct_floor_route_refuted"])
        self.assertTrue(g["strict_task3_closure_not_claimed"])
        self.assertTrue(g["no_qw2191_discharge_claimed"])
        self.assertTrue(g["no_selector_closure_claimed"])
        self.assertTrue(g["no_legacy_kernel_role_transfer"])
        self.assertTrue(g["no_toe_closure_claimed"])


if __name__ == "__main__":
    unittest.main()
