import json
import subprocess
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "fundamental_action_reconstruction" / "scratch" / "bridge_strict_alpha_d12_fourier_phase_reference_obstruction_probe.py"
REPORT_JSON = ROOT / "fundamental_action_reconstruction" / "scratch" / "bridge_strict_alpha_d12_fourier_phase_reference_obstruction_report.json"
REPORT_MD = ROOT / "fundamental_action_reconstruction" / "scratch" / "bridge_strict_alpha_d12_fourier_phase_reference_obstruction_report.md"


class StrictAlphaD12FourierPhaseReferenceObstructionProbeTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        subprocess.run([sys.executable, str(SCRIPT)], check=True, cwd=ROOT)
        cls.payload = json.loads(REPORT_JSON.read_text(encoding="utf-8"))

    def test_result_kind_and_target_replay(self):
        payload = self.payload
        self.assertEqual(
            payload["result_kind"],
            "SCRATCH_STRICT_ALPHA_D12_FOURIER_PHASE_REFERENCE_OBSTRUCTION_PROBE__NOT_A_THEOREM",
        )
        self.assertEqual(payload["target_identity_replay"]["q_power_product"], "256/243")
        self.assertAlmostEqual(payload["target_identity_replay"]["eta_residual_vs_9_5"], 0.0, places=14)
        self.assertEqual(payload["target_identity_replay"]["forward_assignment"], [2, 2, 2, 1, 1])

    def test_previous_no_go_replay(self):
        replay = self.payload["previous_no_go_replay"]
        self.assertEqual(replay["computed_orbit_size"], 24)
        self.assertEqual(replay["computed_stabilizer_size"], 1)
        self.assertEqual(replay["invariant_feature_packets"], 1)

    def test_fourier_power_constant_but_phase_varies(self):
        scan = self.payload["fourier_orbit_scan"]
        self.assertEqual(scan["row_count"], 24)
        self.assertTrue(scan["all_power_spectra_constant_on_orbit"])
        self.assertTrue(all(count == 1 for count in scan["power_unique_counts_by_mode"].values()))
        self.assertGreater(scan["phase_unique_counts_by_mode"]["1"], 1)
        self.assertGreater(scan["phase_unique_counts_by_mode"]["5"], 1)
        self.assertIn("1", scan["nontrivial_phase_modes"])
        self.assertIn("5", scan["nontrivial_phase_modes"])
        self.assertEqual(scan["representative_power_spectrum_modes_0_to_11"]["0"], 64.0)

    def test_dft_covariance_and_guardrails(self):
        payload = self.payload
        covariance = payload["d12_covariance_law_audit"]
        self.assertEqual(covariance["checked_cases"], 288)
        self.assertLess(covariance["max_error"], 1e-10)
        obstruction = payload["phase_reference_obstruction"]
        self.assertIn("D12-invariant", obstruction["invariant_part"])
        self.assertIn("D12-covariant", obstruction["covariant_part"])
        self.assertIn("extra selector data", obstruction["selector_consequence"])
        self.assertEqual(
            payload["candidate_interpretation"]["status"],
            "candidate-supported-but-phase-reference-not-derived",
        )
        hard_limits = "\n".join(payload["hard_limits"])
        self.assertIn("No identity K_legacy_ont == K_strict_gate", hard_limits)
        self.assertIn("D12-covariant", hard_limits)
        self.assertIn("No QW-2191 discharge", hard_limits)
        self.assertIn("No ToE closure", hard_limits)
        self.assertTrue(REPORT_MD.exists())


if __name__ == "__main__":
    unittest.main()
