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
    blob = json.dumps(payload, sort_keys=True, separators=(",", ":"), ensure_ascii=False)
    return hashlib.sha256(blob.encode("utf-8")).hexdigest()


class TestP2318S1268StrictSelectorLaneToTask3MarginBridgeObligationProbe(unittest.TestCase):
    def test_packet(self) -> None:
        subprocess.run(
            [sys.executable, str(ROOT / "p2318_s1268_strict_selector_lane_to_task3_margin_bridge_obligation_probe.py")],
            check=True,
        )
        data = json.loads(
            (G / "p2318_s1268_strict_selector_lane_to_task3_margin_bridge_obligation_probe.json").read_text(
                encoding="utf-8"
            )
        )
        self.assertEqual(data["schema_version"], "p2318_s1268_v1")
        self.assertEqual(data["packet_id"], "P2318")
        self.assertEqual(data["result_kind"], "STRICT_SELECTOR_LANE_TO_PROVIDER_LIFT_OBLIGATION_AUDIT_NO_G1_G3_UPDATE")
        probe = data["strict_selector_lane_to_task3_margin_bridge_obligation_probe"]
        self.assertGreaterEqual(probe["far_bridge_grep_audit"]["hit_count"], 10)
        matrix = probe["interface_obligation_matrix"]
        self.assertEqual(len(matrix), 6)
        self.assertTrue(all(row["status"].startswith("MISSING") for row in matrix))
        diagonal = probe["diagonal_local_lane_scale_witness"]
        shannon = probe["shannon_element_order_lane_scale_witness"]
        self.assertAlmostEqual(diagonal["required_lift"], 0.0068)
        self.assertAlmostEqual(shannon["required_lift"], 0.0068)
        self.assertGreater(diagonal["max_strength"], diagonal["required_lift"])
        self.assertGreater(shannon["max_strength"], shannon["required_lift"])
        self.assertEqual(len(diagonal["witnesses"]), 3)
        self.assertEqual(len(shannon["witnesses"]), 3)
        self.assertTrue(diagonal["witnesses"][0]["passes_required_lift"])
        self.assertFalse(diagonal["witnesses"][1]["passes_required_lift"])
        self.assertFalse(diagonal["witnesses"][2]["passes_required_lift"])
        verdict = probe["bridge_obligation_verdict"]
        self.assertFalse(verdict["target_calibration_is_strict_bridge"])
        self.assertFalse(verdict["all_required_bridge_fields_exported"])
        self.assertEqual(len(verdict["missing_required_bridge_fields"]), 6)
        blockers = probe["existing_bridge_blockers"]
        self.assertTrue(blockers["p2305_norm_to_margin_verdict"]["norm_to_margin_bridge_refuted_as_norm_only_theorem"])
        theorem = probe["theorem_export"]
        self.assertEqual(sha256_json(theorem), probe["theorem_fingerprint_sha256"])
        self.assertIn("do not supply a strict Task-3", theorem["claim"])
        g = data["gatekeeper_checks"]
        self.assertTrue(g["far_bridge_grep_hits_found"])
        self.assertTrue(g["p2302_required_lift_loaded"])
        self.assertTrue(g["p2317_loaded"])
        self.assertTrue(g["diagonal_lane_target_calibration_possible"])
        self.assertTrue(g["shannon_lane_target_calibration_possible"])
        self.assertTrue(g["target_calibration_not_promoted"])
        self.assertTrue(g["missing_bridge_fields_detected"])
        self.assertTrue(g["all_required_bridge_fields_not_exported"])
        self.assertTrue(g["p2305_norm_bridge_refutation_loaded"])
        self.assertTrue(g["strict_g1_g3_not_updated"])
        self.assertTrue(g["no_selector_closure_claimed"])
        self.assertTrue(g["no_qw2191_global_discharge_claimed"])
        self.assertTrue(g["no_toe_closure_claimed"])


if __name__ == "__main__":
    unittest.main()
