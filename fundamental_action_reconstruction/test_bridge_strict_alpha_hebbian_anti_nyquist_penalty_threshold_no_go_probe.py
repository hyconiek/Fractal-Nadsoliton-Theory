import json
import subprocess
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "fundamental_action_reconstruction" / "scratch" / "bridge_strict_alpha_hebbian_anti_nyquist_penalty_threshold_no_go_probe.py"
REPORT_JSON = ROOT / "fundamental_action_reconstruction" / "scratch" / "bridge_strict_alpha_hebbian_anti_nyquist_penalty_threshold_no_go_report.json"
REPORT_MD = ROOT / "fundamental_action_reconstruction" / "scratch" / "bridge_strict_alpha_hebbian_anti_nyquist_penalty_threshold_no_go_report.md"


class StrictAlphaHebbianAntiNyquistPenaltyThresholdNoGoProbeTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        subprocess.run([sys.executable, str(SCRIPT)], check=True, cwd=ROOT)
        cls.payload = json.loads(REPORT_JSON.read_text(encoding="utf-8"))

    def test_result_kind_and_model(self):
        payload = self.payload
        self.assertEqual(
            payload["result_kind"],
            "SCRATCH_STRICT_ALPHA_HEBBIAN_ANTI_NYQUIST_PENALTY_THRESHOLD_NO_GO_PROBE__NOT_A_THEOREM",
        )
        self.assertEqual(payload["status"], "anti-nyquist-penalty-removes-nyquist-but-leaves-k1-k5-unit-tie")
        model = payload["finite_model"]
        self.assertEqual(model["ring"], "Z_12")
        self.assertEqual(model["active_count"], 5)
        self.assertEqual(model["support_count"], 792)
        self.assertEqual(model["histogram_class_count"], 35)
        self.assertEqual(model["mode_histogram_candidate_count"], 210)
        self.assertEqual(model["target_q_power"], "256/243")
        self.assertEqual(model["target_eta"], "9/5")

    def test_threshold_certificate(self):
        threshold = self.payload["threshold_certificate"]
        self.assertEqual(threshold["critical_mu_exact"], {"rational_part": "3/4", "sqrt3_coefficient": "-1/6", "expression": "3/4 - 1/6*sqrt3"})
        self.assertEqual(threshold["d5_reference"]["power_exact"], {"rational_part": "7", "sqrt3_coefficient": "4", "expression": "7 + 4*sqrt3"})
        self.assertEqual(threshold["d5_reference"]["nyquist_power_exact"], "1")
        self.assertEqual(len(threshold["critical_mu_source_rows"]), 1)
        self.assertEqual(threshold["critical_mu_source_rows"][0]["mode"], 6)
        self.assertEqual(threshold["critical_mu_source_rows"][0]["distance_histogram_d1_to_d6"], [0, 4, 0, 4, 0, 2])
        self.assertEqual(threshold["critical_mu_source_rows"][0]["power_exact"], {"rational_part": "25", "sqrt3_coefficient": "0", "expression": "25"})
        self.assertTrue(threshold["all_mu_tie_is_k1_contiguous"])
        self.assertEqual(threshold["all_mu_tie_rows_with_d5"][0]["mode"], 1)
        self.assertEqual(threshold["all_mu_tie_rows_with_d5"][0]["distance_histogram_d1_to_d6"], [4, 3, 2, 1, 0, 0])

    def test_sample_mu_winners(self):
        samples = self.payload["sample_mu_winners"]
        self.assertEqual(samples["mu_0"]["winner_count_over_histogram_mode_candidates"], 1)
        self.assertTrue(samples["mu_0"]["winners"][0]["is_nyquist_sixth"])
        self.assertEqual(samples["mu_critical"]["winner_count_over_histogram_mode_candidates"], 3)
        self.assertEqual(
            sorted((winner["mode"], winner["distance_histogram_d1_to_d6"]) for winner in samples["mu_critical"]["winners"]),
            [
                (1, [4, 3, 2, 1, 0, 0]),
                (5, [0, 3, 2, 1, 4, 0]),
                (6, [0, 4, 0, 4, 0, 2]),
            ],
        )
        self.assertEqual(samples["mu_one_half"]["winner_count_over_histogram_mode_candidates"], 2)
        self.assertTrue(any(winner["is_d5_fifth"] for winner in samples["mu_one_half"]["winners"]))
        self.assertTrue(any(winner["is_contiguous_first"] for winner in samples["mu_one_half"]["winners"]))
        self.assertEqual(samples["mu_1"]["winner_count_over_histogram_mode_candidates"], 2)

    def test_parity_and_guardrails(self):
        parity = self.payload["parity_balance_certificate"]
        self.assertEqual(parity["p6_support_histogram"], {"1": 600, "9": 180, "25": 12})
        self.assertEqual(parity["minimal_p6_exact"], "1")
        self.assertEqual(parity["minimal_p6_support_count"], 600)
        self.assertEqual(parity["minimal_p6_histogram_class_count"], 25)
        self.assertTrue(parity["d5_has_minimal_p6"])
        self.assertTrue(parity["contiguous_has_minimal_p6"])

        readout = self.payload["selector_readout"]
        self.assertIn("k=6 Nyquist winner is removed", readout["positive_filter"])
        self.assertIn("k=1 contiguous", readout["remaining_obstruction"])
        self.assertIn("not derived here", readout["conditional_completion"])

        ontology = self.payload["ontology_guardrail"]
        self.assertIn("nadsoliton itself", ontology["allowed_reading"])
        self.assertIn("No separate informational layer", ontology["forbidden_reading"])

        hard_limits = "\n".join(self.payload["hard_limits"])
        self.assertIn("No identity K_legacy_ont == K_strict_gate", hard_limits)
        self.assertIn("does not uniquely select d5", hard_limits)
        self.assertIn("No QW-2191 discharge", hard_limits)
        self.assertIn("No ToE closure", hard_limits)
        self.assertTrue(REPORT_MD.exists())


if __name__ == "__main__":
    unittest.main()
