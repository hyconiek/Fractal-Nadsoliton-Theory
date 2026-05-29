import json
import subprocess
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "fundamental_action_reconstruction" / "scratch" / "bridge_strict_alpha_self_recorded_arrow_action_lexicographic_selector_probe.py"
REPORT_JSON = ROOT / "fundamental_action_reconstruction" / "scratch" / "bridge_strict_alpha_self_recorded_arrow_action_lexicographic_selector_report.json"
REPORT_MD = ROOT / "fundamental_action_reconstruction" / "scratch" / "bridge_strict_alpha_self_recorded_arrow_action_lexicographic_selector_report.md"


class StrictAlphaSelfRecordedArrowActionLexicographicSelectorProbeTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        subprocess.run([sys.executable, str(SCRIPT)], check=True, cwd=ROOT)
        cls.payload = json.loads(REPORT_JSON.read_text(encoding="utf-8"))

    def test_result_kind_and_target_replay(self):
        payload = self.payload
        self.assertEqual(
            payload["result_kind"],
            "SCRATCH_STRICT_ALPHA_SELF_RECORDED_ARROW_ACTION_LEXICOGRAPHIC_SELECTOR_PROBE__NOT_A_THEOREM",
        )
        self.assertEqual(payload["target_identity_replay"]["q_power_product"], "256/243")
        self.assertAlmostEqual(payload["target_identity_replay"]["eta_residual_vs_9_5"], 0.0, places=14)
        self.assertEqual(payload["target_identity_replay"]["balanced_ledger"], [2, 2, 2, 1, 1])

    def test_previous_stack_replay_and_action_definition(self):
        payload = self.payload
        replay = payload["previous_monotone_stack_replay"]
        self.assertEqual(replay["self_recorded_monotone_count"], 3)
        self.assertTrue(replay["all_checked_power_costs_uniquely_select_balanced"])
        action = payload["action_definition"]
        self.assertIn("sum_i e_i^2", action["ripple_term_R"])
        self.assertIn("max(0, e_{i+1}-e_i)", action["arrow_term_A"])
        self.assertIn("First minimize R", action["lexicographic_action"])

    def test_finite_action_scan(self):
        scan = self.payload["finite_action_scan"]
        self.assertEqual(scan["total_positive_ordered_ledgers"], 35)
        self.assertEqual(scan["ripple_only_minimizers"]["minimum"], 14)
        self.assertEqual(scan["ripple_only_minimizers"]["winner_count"], 10)
        self.assertEqual(scan["arrow_only_minimizers"]["minimum"], 0)
        self.assertEqual(scan["arrow_only_minimizers"]["winner_count"], 3)
        self.assertEqual(scan["lexicographic_ripple_then_arrow_minimizers"]["minimum_pair"], [14, 0])
        self.assertEqual(scan["lexicographic_ripple_then_arrow_minimizers"]["winners"], [[2, 2, 2, 1, 1]])
        self.assertTrue(scan["lexicographic_uniquely_selects_balanced"])
        self.assertEqual(scan["parseval_minimizers"]["winner_count"], 10)

    def test_dominance_certificate_and_guardrails(self):
        payload = self.payload
        dominance = payload["finite_action_scan"]["dominance_certificate"]
        self.assertEqual(dominance["target_ledger"], [2, 2, 2, 1, 1])
        self.assertEqual(dominance["target_ripple"], 14)
        self.assertEqual(dominance["target_arrow_penalty"], 0)
        self.assertEqual(dominance["equal_ripple_competitor_count"], 9)
        self.assertEqual(dominance["higher_ripple_competitor_count"], 25)
        self.assertTrue(dominance["all_equal_ripple_competitors_have_positive_arrow_penalty"])
        self.assertTrue(dominance["all_other_ledgers_have_higher_ripple_or_positive_arrow_tiebreak"])
        self.assertEqual(
            payload["candidate_interpretation"]["status"],
            "candidate-supported-but-arrow-action-not-derived",
        )
        hard_limits = "\n".join(payload["hard_limits"])
        self.assertIn("No identity K_legacy_ont == K_strict_gate", hard_limits)
        self.assertIn("nadsoliton is treated as information", hard_limits)
        self.assertIn("arrow-increase penalty", hard_limits)
        self.assertIn("No QW-2191 discharge", hard_limits)
        self.assertIn("No ToE closure", hard_limits)
        self.assertTrue(REPORT_MD.exists())


if __name__ == "__main__":
    unittest.main()
