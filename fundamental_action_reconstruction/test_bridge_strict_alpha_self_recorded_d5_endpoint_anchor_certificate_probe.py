import json
import subprocess
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "fundamental_action_reconstruction" / "scratch" / "bridge_strict_alpha_self_recorded_d5_endpoint_anchor_certificate_probe.py"
REPORT_JSON = ROOT / "fundamental_action_reconstruction" / "scratch" / "bridge_strict_alpha_self_recorded_d5_endpoint_anchor_certificate_report.json"
REPORT_MD = ROOT / "fundamental_action_reconstruction" / "scratch" / "bridge_strict_alpha_self_recorded_d5_endpoint_anchor_certificate_report.md"


class StrictAlphaSelfRecordedD5EndpointAnchorCertificateProbeTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        subprocess.run([sys.executable, str(SCRIPT)], check=True, cwd=ROOT)
        cls.payload = json.loads(REPORT_JSON.read_text(encoding="utf-8"))

    def test_result_kind_and_target_replay(self):
        payload = self.payload
        self.assertEqual(
            payload["result_kind"],
            "SCRATCH_STRICT_ALPHA_SELF_RECORDED_D5_ENDPOINT_ANCHOR_CERTIFICATE_PROBE__NOT_A_THEOREM",
        )
        self.assertEqual(payload["target_identity_replay"]["q_power_product"], "256/243")
        self.assertAlmostEqual(payload["target_identity_replay"]["eta_residual_vs_9_5"], 0.0, places=14)
        self.assertEqual(payload["target_identity_replay"]["forward_assignment"], [2, 2, 2, 1, 1])

    def test_previous_phase_aliasing_replay_and_rule(self):
        payload = self.payload
        replay = payload["previous_phase_aliasing_replay"]
        self.assertEqual(replay["source_complete_modes"], [1, 5, 7, 11])
        self.assertTrue(replay["all_observed_counts_match_gcd_formula"])
        rule = payload["self_recorded_anchor_rule"]
        self.assertIn("unique value-2 endpoint", rule["rule"])
        self.assertIn("internal endpoint record", rule["ontology_read"])
        self.assertIn("D12-equivariant", rule["not_an_absolute_origin"])

    def test_anchor_reconstruction_scan(self):
        scan = self.payload["anchor_reconstruction_scan"]
        self.assertEqual(scan["row_count"], 24)
        self.assertTrue(scan["all_sources_recovered"])
        self.assertTrue(scan["all_orientations_recovered"])
        self.assertTrue(scan["all_ordered_values_match_forward_assignment"])
        sample = next(
            row for row in scan["rows"]
            if row["actual_source"] == 7 and row["actual_orientation"] == -1
        )
        self.assertEqual(sample["inferred_anchor"]["inferred_source"], 7)
        self.assertEqual(sample["inferred_anchor"]["inferred_orientation"], -1)
        self.assertEqual(sample["inferred_anchor"]["ordered_values"], [2, 2, 2, 1, 1])
        self.assertEqual(sorted(sample["inferred_anchor"]["endpoint_values"].values()), [1, 2])

    def test_equivariance_and_guardrails(self):
        payload = self.payload
        equivariance = payload["d12_equivariance_audit"]
        self.assertEqual(equivariance["checked_cases"], 576)
        self.assertEqual(equivariance["mismatch_count"], 0)
        self.assertTrue(equivariance["all_cases_equivariant"])
        consequence = payload["selector_consequence"]
        self.assertIn("self-recorded", consequence["what_is_gained"])
        self.assertIn("equivariant", consequence["why_this_does_not_contradict_D12_no_go"])
        self.assertIn("d5 support", consequence["what_remains_unproved"])
        self.assertEqual(
            payload["candidate_interpretation"]["status"],
            "candidate-supported-but-d5-support-and-ledger-not-derived",
        )
        hard_limits = "\n".join(payload["hard_limits"])
        self.assertIn("No identity K_legacy_ont == K_strict_gate", hard_limits)
        self.assertIn("nadsoliton is treated as information", hard_limits)
        self.assertIn("D12-equivariant self-record extractor", hard_limits)
        self.assertIn("No QW-2191 discharge", hard_limits)
        self.assertIn("No ToE closure", hard_limits)
        self.assertTrue(REPORT_MD.exists())


if __name__ == "__main__":
    unittest.main()
