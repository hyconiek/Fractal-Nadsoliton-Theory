import json
import subprocess
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "fundamental_action_reconstruction" / "scratch" / "bridge_strict_alpha_self_recorded_monotone_min_ripple_stack_probe.py"
REPORT_JSON = ROOT / "fundamental_action_reconstruction" / "scratch" / "bridge_strict_alpha_self_recorded_monotone_min_ripple_stack_report.json"
REPORT_MD = ROOT / "fundamental_action_reconstruction" / "scratch" / "bridge_strict_alpha_self_recorded_monotone_min_ripple_stack_report.md"


class StrictAlphaSelfRecordedMonotoneMinRippleStackProbeTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        subprocess.run([sys.executable, str(SCRIPT)], check=True, cwd=ROOT)
        cls.payload = json.loads(REPORT_JSON.read_text(encoding="utf-8"))

    def test_result_kind_and_target_replay(self):
        payload = self.payload
        self.assertEqual(
            payload["result_kind"],
            "SCRATCH_STRICT_ALPHA_SELF_RECORDED_MONOTONE_MIN_RIPPLE_STACK_PROBE__NOT_A_THEOREM",
        )
        self.assertEqual(payload["target_identity_replay"]["q_power_product"], "256/243")
        self.assertAlmostEqual(payload["target_identity_replay"]["eta_residual_vs_9_5"], 0.0, places=14)
        self.assertEqual(payload["target_identity_replay"]["balanced_ledger"], [2, 2, 2, 1, 1])

    def test_previous_nonuniqueness_replay_and_stack_definition(self):
        payload = self.payload
        replay = payload["previous_nonuniqueness_replay"]
        self.assertEqual(replay["total_positive_ordered_ledgers"], 35)
        self.assertEqual(replay["endpoint_distinct_self_recording_count"], 22)
        definition = payload["selector_stack_definition"]
        self.assertIn("e_0 > e_4", definition["stage_1_self_record"])
        self.assertIn("monotone non-increasing", definition["stage_2_monotone_formation"])
        self.assertIn("selector premises", definition["premise_warning"])

    def test_stack_scan_and_minimizers(self):
        scan = self.payload["finite_stack_scan"]
        self.assertEqual(scan["total_positive_ordered_ledgers"], 35)
        self.assertEqual(scan["endpoint_self_recorded_from_left_count"], 11)
        self.assertEqual(scan["self_recorded_monotone_count"], 3)
        self.assertEqual(
            scan["self_recorded_monotone_ledgers"],
            [[2, 2, 2, 1, 1], [3, 2, 1, 1, 1], [4, 1, 1, 1, 1]],
        )
        self.assertEqual(scan["sum_square_minimizers_on_monotone_stack"]["winners"], [[2, 2, 2, 1, 1]])
        self.assertEqual(scan["parseval_minimizers_on_monotone_stack"]["winners"], [[2, 2, 2, 1, 1]])
        self.assertTrue(scan["all_checked_power_costs_uniquely_select_balanced"])
        for packet in scan["power_minimizers_on_monotone_stack"].values():
            self.assertEqual(packet["winner_count"], 1)
            self.assertEqual(packet["winners"], [[2, 2, 2, 1, 1]])

    def test_smoothing_certificate_and_guardrails(self):
        payload = self.payload
        smoothing = payload["smoothing_certificate_on_monotone_stack"]
        self.assertEqual(len(smoothing), 2)
        self.assertEqual(smoothing[0]["sum_squares_drop"], 4)
        self.assertEqual(smoothing[1]["sum_squares_drop"], 2)
        self.assertEqual(smoothing[0]["parseval_non_dc_drop"], 20)
        self.assertEqual(smoothing[1]["parseval_non_dc_drop"], 10)
        consequence = payload["selector_consequence"]
        self.assertIn("three canonical descending ledgers", consequence["what_self_record_plus_monotone_does"])
        self.assertIn("uniquely selects", consequence["what_min_ripple_adds"])
        self.assertEqual(
            payload["candidate_interpretation"]["status"],
            "candidate-supported-but-monotone-min-ripple-premises-not-derived",
        )
        hard_limits = "\n".join(payload["hard_limits"])
        self.assertIn("No identity K_legacy_ont == K_strict_gate", hard_limits)
        self.assertIn("nadsoliton is treated as information", hard_limits)
        self.assertIn("conditional selector stack", hard_limits)
        self.assertIn("No QW-2191 discharge", hard_limits)
        self.assertIn("No ToE closure", hard_limits)
        self.assertTrue(REPORT_MD.exists())


if __name__ == "__main__":
    unittest.main()
