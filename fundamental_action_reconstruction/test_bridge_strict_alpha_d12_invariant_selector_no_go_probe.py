import json
import subprocess
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "fundamental_action_reconstruction" / "scratch" / "bridge_strict_alpha_d12_invariant_selector_no_go_probe.py"
REPORT_JSON = ROOT / "fundamental_action_reconstruction" / "scratch" / "bridge_strict_alpha_d12_invariant_selector_no_go_report.json"
REPORT_MD = ROOT / "fundamental_action_reconstruction" / "scratch" / "bridge_strict_alpha_d12_invariant_selector_no_go_report.md"


class StrictAlphaD12InvariantSelectorNoGoProbeTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        subprocess.run([sys.executable, str(SCRIPT)], check=True, cwd=ROOT)
        cls.payload = json.loads(REPORT_JSON.read_text(encoding="utf-8"))

    def test_result_kind_and_target_replay(self):
        payload = self.payload
        self.assertEqual(
            payload["result_kind"],
            "SCRATCH_STRICT_ALPHA_D12_INVARIANT_SELECTOR_NO_GO_PROBE__NOT_A_THEOREM",
        )
        self.assertEqual(payload["target_identity_replay"]["q_power_product"], "256/243")
        self.assertAlmostEqual(payload["target_identity_replay"]["eta_residual_vs_9_5"], 0.0, places=14)
        self.assertEqual(payload["target_identity_replay"]["forward_assignment"], [2, 2, 2, 1, 1])

    def test_previous_orbit_replay_and_no_go_statement(self):
        payload = self.payload
        replay = payload["previous_symmetry_audit_replay"]
        self.assertTrue(replay["anchored_rows_equal_full_dihedral_orbit"])
        self.assertTrue(replay["value_configuration_has_trivial_stabilizer"])
        self.assertEqual(replay["reported_value_configuration_orbit_size"], 24)
        self.assertEqual(replay["reported_unoriented_support_orbit_size"], 12)
        orbit = payload["d12_orbit_replay"]
        self.assertEqual(orbit["dihedral_group_order"], 24)
        self.assertEqual(orbit["computed_orbit_size"], 24)
        self.assertEqual(orbit["computed_stabilizer_size"], 1)
        self.assertIn("constant", payload["finite_no_go_statement"]["conclusion"])
        self.assertIn("non-invariant", payload["finite_no_go_statement"]["what_can_select"])

    def test_invariant_feature_collapse(self):
        audit = self.payload["invariant_feature_audit"]
        self.assertEqual(audit["row_count"], 24)
        self.assertEqual(audit["unique_invariant_feature_packet_count"], 1)
        self.assertTrue(audit["all_checked_invariant_features_constant"])
        self.assertTrue(all(count == 1 for count in audit["feature_class_unique_counts"].values()))
        representative = audit["representative_invariant_feature_packet"]
        self.assertEqual(representative["value_count_vector"], [7, 2, 3])
        self.assertEqual(representative["sorted_nonzero_values_desc"], [2, 2, 2, 1, 1])
        self.assertEqual(representative["support_distance_histogram_d1_to_d6"], [0, 3, 2, 1, 4, 0])
        self.assertEqual(representative["total_mass"], 8)
        self.assertEqual(representative["total_square_mass"], 14)

    def test_noninvariant_contrast_and_guardrails(self):
        payload = self.payload
        fixed = payload["noninvariant_contrast"]["fixed_frame_score_summary"]
        self.assertGreater(fixed["linear_node_moment"]["unique_value_count"], 1)
        self.assertGreater(fixed["quadratic_node_moment"]["unique_value_count"], 1)
        self.assertGreater(fixed["right_half_mass_nodes_6_to_11"]["unique_value_count"], 1)
        self.assertEqual(
            payload["candidate_interpretation"]["status"],
            "candidate-supported-but-D12-invariant-selectors-provably-insufficient",
        )
        hard_limits = "\n".join(payload["hard_limits"])
        self.assertIn("No identity K_legacy_ont == K_strict_gate", hard_limits)
        self.assertIn("No D12-invariant selector", hard_limits)
        self.assertIn("No QW-2191 discharge", hard_limits)
        self.assertIn("No ToE closure", hard_limits)
        self.assertTrue(REPORT_MD.exists())


if __name__ == "__main__":
    unittest.main()
