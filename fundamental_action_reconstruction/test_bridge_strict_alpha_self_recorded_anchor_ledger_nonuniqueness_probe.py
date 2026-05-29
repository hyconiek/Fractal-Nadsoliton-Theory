import json
import subprocess
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "fundamental_action_reconstruction" / "scratch" / "bridge_strict_alpha_self_recorded_anchor_ledger_nonuniqueness_probe.py"
REPORT_JSON = ROOT / "fundamental_action_reconstruction" / "scratch" / "bridge_strict_alpha_self_recorded_anchor_ledger_nonuniqueness_report.json"
REPORT_MD = ROOT / "fundamental_action_reconstruction" / "scratch" / "bridge_strict_alpha_self_recorded_anchor_ledger_nonuniqueness_report.md"


class StrictAlphaSelfRecordedAnchorLedgerNonuniquenessProbeTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        subprocess.run([sys.executable, str(SCRIPT)], check=True, cwd=ROOT)
        cls.payload = json.loads(REPORT_JSON.read_text(encoding="utf-8"))

    def test_result_kind_and_target_replay(self):
        payload = self.payload
        self.assertEqual(
            payload["result_kind"],
            "SCRATCH_STRICT_ALPHA_SELF_RECORDED_ANCHOR_LEDGER_NONUNIQUENESS_PROBE__NOT_A_THEOREM",
        )
        self.assertEqual(payload["target_identity_replay"]["q_power_product"], "256/243")
        self.assertAlmostEqual(payload["target_identity_replay"]["eta_residual_vs_9_5"], 0.0, places=14)
        self.assertEqual(payload["target_identity_replay"]["balanced_ledger"], [2, 2, 2, 1, 1])

    def test_previous_self_recorded_replay_and_criterion(self):
        payload = self.payload
        replay = payload["previous_self_recorded_anchor_replay"]
        self.assertTrue(replay["all_sources_recovered"])
        self.assertTrue(replay["all_orientations_recovered"])
        self.assertEqual(replay["d12_equivariance_checked_cases"], 576)
        criterion = payload["endpoint_anchor_criterion"]
        self.assertIn("endpoint values are unequal", criterion["criterion"])
        self.assertIn("Equal endpoint values", criterion["failure_mode"])
        self.assertIn("does not rank", criterion["why_not_selector"])

    def test_ordered_ledger_counts_and_partition_summary(self):
        payload = self.payload
        scan = payload["ordered_ledger_scan"]
        self.assertEqual(scan["total_positive_ordered_ledgers"], 35)
        self.assertEqual(scan["endpoint_distinct_self_recording_count"], 22)
        self.assertEqual(scan["endpoint_equal_ambiguous_count"], 13)
        summary = payload["canonical_partition_summary"]
        self.assertEqual(summary["4,1,1,1,1"]["ordered_permutation_count"], 5)
        self.assertEqual(summary["4,1,1,1,1"]["endpoint_distinct_count"], 2)
        self.assertEqual(summary["3,2,1,1,1"]["endpoint_distinct_count"], 14)
        self.assertEqual(summary["2,2,2,1,1"]["endpoint_distinct_count"], 6)
        self.assertTrue(all(packet["has_endpoint_anchor_examples"] for packet in summary.values()))

    def test_canonical_replay_and_guardrails(self):
        payload = self.payload
        replay = payload["canonical_partition_replay"]
        for key in ("4,1,1,1,1", "3,2,1,1,1", "2,2,2,1,1"):
            self.assertTrue(replay[key]["sorted_descending_order_endpoint_anchor_status"]["endpoint_values_distinct"])
        consequence = payload["selector_consequence"]
        self.assertIn("not select the balanced ledger", consequence["what_is_ruled_out"])
        self.assertIn("ledger selector", consequence["remaining_selector_obligation"])
        self.assertEqual(
            payload["candidate_interpretation"]["status"],
            "candidate-supported-but-self-record-anchor-not-ledger-selector",
        )
        hard_limits = "\n".join(payload["hard_limits"])
        self.assertIn("No identity K_legacy_ont == K_strict_gate", hard_limits)
        self.assertIn("nadsoliton is treated as information", hard_limits)
        self.assertIn("not a ledger selector theorem", hard_limits)
        self.assertIn("No QW-2191 discharge", hard_limits)
        self.assertIn("No ToE closure", hard_limits)
        self.assertTrue(REPORT_MD.exists())


if __name__ == "__main__":
    unittest.main()
