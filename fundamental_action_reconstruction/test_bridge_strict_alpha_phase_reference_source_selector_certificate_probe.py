import json
import subprocess
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "fundamental_action_reconstruction" / "scratch" / "bridge_strict_alpha_phase_reference_source_selector_certificate_probe.py"
REPORT_JSON = ROOT / "fundamental_action_reconstruction" / "scratch" / "bridge_strict_alpha_phase_reference_source_selector_certificate_report.json"
REPORT_MD = ROOT / "fundamental_action_reconstruction" / "scratch" / "bridge_strict_alpha_phase_reference_source_selector_certificate_report.md"


class StrictAlphaPhaseReferenceSourceSelectorCertificateProbeTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        subprocess.run([sys.executable, str(SCRIPT)], check=True, cwd=ROOT)
        cls.payload = json.loads(REPORT_JSON.read_text(encoding="utf-8"))

    def test_result_kind_and_target_replay(self):
        payload = self.payload
        self.assertEqual(
            payload["result_kind"],
            "SCRATCH_STRICT_ALPHA_PHASE_REFERENCE_SOURCE_SELECTOR_CERTIFICATE_PROBE__NOT_A_THEOREM",
        )
        self.assertEqual(payload["target_identity_replay"]["q_power_product"], "256/243")
        self.assertAlmostEqual(payload["target_identity_replay"]["eta_residual_vs_9_5"], 0.0, places=14)
        self.assertEqual(payload["target_identity_replay"]["forward_assignment"], [2, 2, 2, 1, 1])

    def test_previous_no_go_replay_and_formula(self):
        payload = self.payload
        replay = payload["previous_source_no_go_replay"]
        self.assertEqual(replay["source_orbit_minus_size"], 12)
        self.assertEqual(replay["source_orbit_plus_size"], 12)
        self.assertEqual(replay["full_bispectrum_minus_signature_count"], 1)
        self.assertEqual(replay["full_bispectrum_plus_signature_count"], 1)
        formula = payload["phase_recovery_formula"]
        self.assertEqual(formula["required_modes"], [1, 5, 7, 11])
        self.assertIn("inv(k mod 12)", formula["formula"])
        self.assertIn("invertible modulo 12", formula["why_coprime"])
        self.assertEqual(set(formula["calibrated_source_0_phase_references"].keys()), {"-1", "1"})

    def test_conditional_selector_recovers_all_rows(self):
        scan = self.payload["conditional_selector_scan"]
        self.assertEqual(scan["row_count"], 24)
        self.assertTrue(scan["all_orientations_recovered_by_chiral_bispectrum"])
        self.assertTrue(scan["all_rows_recovered_by_all_coprime_modes"])
        self.assertTrue(all(scan["source_recovery_success_by_mode"].values()))
        self.assertEqual(scan["unique_predicted_orientation_source_pairs_using_mode_1"], 24)
        self.assertEqual(scan["expected_pair_count"], 24)
        sample = next(
            row for row in scan["rows"]
            if row["actual_source"] == 7 and row["actual_orientation"] == -1
        )
        self.assertEqual(sample["predicted_orientation_from_chiral_bispectrum"], -1)
        self.assertTrue(sample["all_modes_recover_source"])
        for mode in ("1", "5", "7", "11"):
            self.assertEqual(sample["per_coprime_mode_source_recovery"][mode]["recovered_source"], 7)

    def test_premise_based_guardrails(self):
        payload = self.payload
        consequence = payload["selector_consequence"]
        self.assertIn("recovers orientation and source", consequence["what_is_sufficient"])
        self.assertIn("external premises", consequence["what_is_not_derived"])
        self.assertIn("does not contradict", consequence["relation_to_no_go"])
        self.assertEqual(
            payload["candidate_interpretation"]["status"],
            "candidate-supported-but-phase-origin-premise-not-derived",
        )
        hard_limits = "\n".join(payload["hard_limits"])
        self.assertIn("No identity K_legacy_ont == K_strict_gate", hard_limits)
        self.assertIn("premise", hard_limits)
        self.assertIn("phase-origin", hard_limits)
        self.assertIn("No QW-2191 discharge", hard_limits)
        self.assertIn("No ToE closure", hard_limits)
        self.assertTrue(REPORT_MD.exists())


if __name__ == "__main__":
    unittest.main()
