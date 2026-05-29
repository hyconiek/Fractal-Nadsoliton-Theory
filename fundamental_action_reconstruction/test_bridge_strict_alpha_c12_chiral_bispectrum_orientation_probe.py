import json
import subprocess
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "fundamental_action_reconstruction" / "scratch" / "bridge_strict_alpha_c12_chiral_bispectrum_orientation_probe.py"
REPORT_JSON = ROOT / "fundamental_action_reconstruction" / "scratch" / "bridge_strict_alpha_c12_chiral_bispectrum_orientation_report.json"
REPORT_MD = ROOT / "fundamental_action_reconstruction" / "scratch" / "bridge_strict_alpha_c12_chiral_bispectrum_orientation_report.md"


class StrictAlphaC12ChiralBispectrumOrientationProbeTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        subprocess.run([sys.executable, str(SCRIPT)], check=True, cwd=ROOT)
        cls.payload = json.loads(REPORT_JSON.read_text(encoding="utf-8"))

    def test_result_kind_and_target_replay(self):
        payload = self.payload
        self.assertEqual(
            payload["result_kind"],
            "SCRATCH_STRICT_ALPHA_C12_CHIRAL_BISPECTRUM_ORIENTATION_PROBE__NOT_A_THEOREM",
        )
        self.assertEqual(payload["target_identity_replay"]["q_power_product"], "256/243")
        self.assertAlmostEqual(payload["target_identity_replay"]["eta_residual_vs_9_5"], 0.0, places=14)
        self.assertEqual(payload["target_identity_replay"]["forward_assignment"], [2, 2, 2, 1, 1])

    def test_previous_phase_obstruction_replay(self):
        replay = self.payload["previous_phase_obstruction_replay"]
        self.assertTrue(replay["all_power_spectra_constant_on_orbit"])
        self.assertLess(replay["dft_covariance_max_error"], 1e-10)
        self.assertEqual(replay["phase_reference_status"], "candidate-supported-but-phase-reference-not-derived")

    def test_bispectrum_orientation_and_source_scan(self):
        scan = self.payload["c12_chiral_orientation_scan"]
        self.assertEqual(scan["row_count"], 24)
        self.assertEqual(scan["orientation_separating_pair_count"], 5)
        self.assertTrue(scan["all_selected_pairs_translation_invariant_over_sources"])
        self.assertTrue(scan["all_selected_pairs_reflection_chiral"])
        pair = scan["orientation_bispectrum_summary"]["1,5"]
        self.assertEqual(pair["by_orientation"]["-1"]["unique_value_count_over_12_sources"], 1)
        self.assertEqual(pair["by_orientation"]["1"]["unique_value_count_over_12_sources"], 1)
        self.assertEqual(pair["by_orientation"]["-1"]["representative_value"], {"real": 12.0, "imag": 2.0})
        self.assertEqual(pair["by_orientation"]["1"]["representative_value"], {"real": 12.0, "imag": -2.0})
        self.assertTrue(pair["orientation_separating"])

    def test_source_degeneracy_and_guardrails(self):
        payload = self.payload
        source = payload["source_degeneracy_after_c12_chiral_phase"]
        self.assertEqual(source["-1"]["signature_count"], 1)
        self.assertEqual(source["1"]["signature_count"], 1)
        self.assertEqual(source["-1"]["source_counts_per_signature"], [12])
        self.assertEqual(source["1"]["source_counts_per_signature"], [12])
        consequence = payload["selector_consequence"]
        self.assertIn("handedness", consequence["what_is_gained"])
        self.assertIn("12-fold", consequence["what_remains_unselected"])
        self.assertIn("flips under reflection", consequence["why_this_is_not_D12_invariant"])
        self.assertEqual(
            payload["candidate_interpretation"]["status"],
            "candidate-supported-but-handedness-and-source-not-derived",
        )
        hard_limits = "\n".join(payload["hard_limits"])
        self.assertIn("No identity K_legacy_ont == K_strict_gate", hard_limits)
        self.assertIn("not D12-invariant", hard_limits)
        self.assertIn("12-fold degenerate", hard_limits)
        self.assertIn("No QW-2191 discharge", hard_limits)
        self.assertIn("No ToE closure", hard_limits)
        self.assertTrue(REPORT_MD.exists())


if __name__ == "__main__":
    unittest.main()
