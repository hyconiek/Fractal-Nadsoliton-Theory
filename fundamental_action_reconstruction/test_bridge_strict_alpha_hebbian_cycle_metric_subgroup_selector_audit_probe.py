import json
import subprocess
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "fundamental_action_reconstruction" / "scratch" / "bridge_strict_alpha_hebbian_cycle_metric_subgroup_selector_audit_probe.py"
REPORT_JSON = ROOT / "fundamental_action_reconstruction" / "scratch" / "bridge_strict_alpha_hebbian_cycle_metric_subgroup_selector_audit_report.json"
REPORT_MD = ROOT / "fundamental_action_reconstruction" / "scratch" / "bridge_strict_alpha_hebbian_cycle_metric_subgroup_selector_audit_report.md"


class StrictAlphaHebbianCycleMetricSubgroupSelectorAuditProbeTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        subprocess.run([sys.executable, str(SCRIPT)], check=True, cwd=ROOT)
        cls.payload = json.loads(REPORT_JSON.read_text(encoding="utf-8"))

    def test_result_kind_and_model(self):
        payload = self.payload
        self.assertEqual(
            payload["result_kind"],
            "SCRATCH_STRICT_ALPHA_HEBBIAN_CYCLE_METRIC_SUBGROUP_SELECTOR_AUDIT_PROBE__NOT_A_THEOREM",
        )
        self.assertEqual(
            payload["status"],
            "cycle-metric-record-supplies-required-subgroup-conditionally-not-strict-source",
        )
        model = payload["finite_model"]
        self.assertEqual(model["ring"], "Z_12")
        self.assertEqual(model["active_count"], 5)
        self.assertEqual(model["support_count"], 792)
        self.assertEqual(model["automorphism_units"], [1, 5, 7, 11])
        self.assertEqual(model["minimal_required_subgroup_from_previous_probe"], [1, 11])
        self.assertEqual(model["target_q_power"], "256/243")
        self.assertEqual(model["target_eta"], "9/5")

    def test_distance_stabilizers(self):
        audit = self.payload["distance_stabilizer_audit"]
        self.assertEqual(
            audit["distance_stabilizers"],
            {"1": [1, 11], "2": [1, 5, 7, 11], "3": [1, 5, 7, 11], "4": [1, 5, 7, 11], "5": [1, 11], "6": [1, 5, 7, 11]},
        )
        self.assertEqual(audit["nearest_neighbor_distance_1_stabilizer"], [1, 11])
        self.assertEqual(audit["d5_step_distance_5_stabilizer"], [1, 11])
        self.assertTrue(audit["distance_1_matches_minimal_required_subgroup"])
        self.assertTrue(audit["distance_5_matches_minimal_required_subgroup"])
        row1 = next(row for row in audit["rows"] if row["folded_distance"] == 1)
        self.assertEqual(row1["unit_edge_preservation"], {"1": True, "5": False, "7": False, "11": True})
        self.assertEqual(row1["edge_count"], 12)

    def test_selector_replay(self):
        replay = self.payload["selector_replay"]
        full_aut = replay["full_aut_selector"]
        reduced = replay["cycle_metric_reduced_selector"]
        self.assertFalse(full_aut["selects_d5_in_reduced_symmetry"])
        self.assertFalse(full_aut["d5_singleton_invariant"])
        self.assertTrue(reduced["selects_d5_in_reduced_symmetry"])
        self.assertTrue(reduced["d5_singleton_invariant"])
        self.assertTrue(reduced["contiguous_singleton_invariant"])
        self.assertEqual(reduced["label_orbits_under_subgroup"], {"contiguous_step_1_11": [1], "fifth_step_d5_step_5_7": [5]})
        self.assertTrue(replay["cycle_metric_makes_d5_selector_invariant"])
        self.assertFalse(replay["full_aut_makes_d5_selector_invariant"])

    def test_interpretation_and_guardrails(self):
        payload = self.payload
        interpretation = payload["candidate_source_interpretation"]
        self.assertIn("exactly the subgroup effect", interpretation["finite_gain"])
        self.assertIn("selects d5 invariantly", interpretation["conditional_positive_result"])
        self.assertIn("not derived", interpretation["honest_limit"])
        self.assertIn("concrete finite candidate", interpretation["relation_to_previous_probe"])

        hard_limits = "\n".join(payload["hard_limits"])
        self.assertIn("No identity K_legacy_ont == K_strict_gate", hard_limits)
        self.assertIn("No theorem derives the d5 teacher/self-record trace", hard_limits)
        self.assertIn("No theorem derives the cycle metric/locality record", hard_limits)
        self.assertIn("conditional subgroup-source audit", hard_limits)
        self.assertIn("No QW-2191 discharge", hard_limits)
        self.assertIn("No ToE closure", hard_limits)
        self.assertTrue(REPORT_MD.exists())


if __name__ == "__main__":
    unittest.main()
