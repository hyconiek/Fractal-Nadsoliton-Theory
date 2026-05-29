import json
import subprocess
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "fundamental_action_reconstruction" / "scratch" / "bridge_strict_alpha_c12_source_translation_no_go_probe.py"
REPORT_JSON = ROOT / "fundamental_action_reconstruction" / "scratch" / "bridge_strict_alpha_c12_source_translation_no_go_report.json"
REPORT_MD = ROOT / "fundamental_action_reconstruction" / "scratch" / "bridge_strict_alpha_c12_source_translation_no_go_report.md"


class StrictAlphaC12SourceTranslationNoGoProbeTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        subprocess.run([sys.executable, str(SCRIPT)], check=True, cwd=ROOT)
        cls.payload = json.loads(REPORT_JSON.read_text(encoding="utf-8"))

    def test_result_kind_and_target_replay(self):
        payload = self.payload
        self.assertEqual(
            payload["result_kind"],
            "SCRATCH_STRICT_ALPHA_C12_SOURCE_TRANSLATION_NO_GO_PROBE__NOT_A_THEOREM",
        )
        self.assertEqual(payload["target_identity_replay"]["q_power_product"], "256/243")
        self.assertAlmostEqual(payload["target_identity_replay"]["eta_residual_vs_9_5"], 0.0, places=14)
        self.assertEqual(payload["target_identity_replay"]["forward_assignment"], [2, 2, 2, 1, 1])

    def test_previous_chiral_replay_and_no_go_statement(self):
        payload = self.payload
        replay = payload["previous_chiral_bispectrum_replay"]
        self.assertEqual(replay["orientation_separating_pair_count"], 5)
        self.assertEqual(replay["source_degeneracy_minus"], [12])
        self.assertEqual(replay["source_degeneracy_plus"], [12])
        statement = payload["finite_no_go_statement"]
        self.assertIn("C12-translation-invariant", statement["premise"])
        self.assertIn("constant", statement["conclusion"])
        self.assertIn("origin/phase reference", statement["what_can_select"])

    def test_source_orbits_and_invariant_signatures(self):
        audits = self.payload["source_orbit_audits_by_orientation"]
        for orientation in ("-1", "1"):
            audit = audits[orientation]
            self.assertEqual(audit["computed_translation_orbit_size"], 12)
            self.assertEqual(audit["computed_translation_stabilizer_size"], 1)
            self.assertEqual(audit["computed_translation_stabilizer"], [0])
            self.assertEqual(audit["row_count"], 12)
            self.assertEqual(audit["power_signature_unique_count_over_sources"], 1)
            self.assertEqual(audit["full_bispectrum_signature_unique_count_over_sources"], 1)
            self.assertGreater(audit["raw_phase_signature_unique_count_over_sources"], 1)
            self.assertGreater(len(audit["raw_phase_mode_1_unique_values"]), 1)
        self.assertEqual(audits["-1"]["representative_bispectrum_sample"]["1,5"], {"real": 12.0, "imag": 2.0})
        self.assertEqual(audits["1"]["representative_bispectrum_sample"]["1,5"], {"real": 12.0, "imag": -2.0})

    def test_translation_covariance_and_guardrails(self):
        payload = self.payload
        covariance = payload["translation_covariance_and_invariance_audit"]
        self.assertEqual(covariance["checked_dft_cases"], 288)
        self.assertEqual(covariance["checked_power_cases"], 288)
        self.assertEqual(covariance["checked_bispectrum_cases"], 3456)
        self.assertLess(covariance["max_dft_covariance_error"], 1e-10)
        self.assertLess(covariance["max_power_invariance_error"], 1e-10)
        self.assertLess(covariance["max_bispectrum_invariance_error"], 1e-10)
        self.assertEqual(
            payload["candidate_interpretation"]["status"],
            "candidate-supported-but-source-localizer-not-derived",
        )
        hard_limits = "\n".join(payload["hard_limits"])
        self.assertIn("No identity K_legacy_ont == K_strict_gate", hard_limits)
        self.assertIn("C12-translation-invariant scores cannot choose", hard_limits)
        self.assertIn("Raw Fourier phase variation requires", hard_limits)
        self.assertIn("No QW-2191 discharge", hard_limits)
        self.assertIn("No ToE closure", hard_limits)
        self.assertTrue(REPORT_MD.exists())


if __name__ == "__main__":
    unittest.main()
