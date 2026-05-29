import json
import subprocess
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "fundamental_action_reconstruction" / "scratch" / "bridge_strict_alpha_hebbian_minimal_aut_breaking_selector_probe.py"
REPORT_JSON = ROOT / "fundamental_action_reconstruction" / "scratch" / "bridge_strict_alpha_hebbian_minimal_aut_breaking_selector_report.json"
REPORT_MD = ROOT / "fundamental_action_reconstruction" / "scratch" / "bridge_strict_alpha_hebbian_minimal_aut_breaking_selector_report.md"


class StrictAlphaHebbianMinimalAutBreakingSelectorProbeTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        subprocess.run([sys.executable, str(SCRIPT)], check=True, cwd=ROOT)
        cls.payload = json.loads(REPORT_JSON.read_text(encoding="utf-8"))

    def test_result_kind_and_model(self):
        payload = self.payload
        self.assertEqual(
            payload["result_kind"],
            "SCRATCH_STRICT_ALPHA_HEBBIAN_MINIMAL_AUT_BREAKING_SELECTOR_PROBE__NOT_A_THEOREM",
        )
        self.assertEqual(payload["status"], "minimal-subgroup-reduction-for-d5-label-selector-not-strict-source")
        model = payload["finite_model"]
        self.assertEqual(model["ring"], "Z_12")
        self.assertEqual(model["active_count"], 5)
        self.assertEqual(model["support_count"], 792)
        self.assertEqual(model["automorphism_units"], [1, 5, 7, 11])
        self.assertEqual(model["target_q_power"], "256/243")
        self.assertEqual(model["target_eta"], "9/5")

    def test_subgroup_enumeration(self):
        enum = self.payload["aut_subgroup_enumeration"]
        self.assertEqual(enum["subgroup_count"], 5)
        self.assertEqual(enum["subgroups"], [[1], [1, 5], [1, 7], [1, 11], [1, 5, 7, 11]])
        rows = {tuple(row["subgroup"]): row for row in enum["rows"]}
        self.assertFalse(rows[(1, 5, 7, 11)]["d5_singleton_invariant"])
        self.assertFalse(rows[(1, 5)]["d5_singleton_invariant"])
        self.assertFalse(rows[(1, 7)]["d5_singleton_invariant"])
        self.assertTrue(rows[(1, 11)]["d5_singleton_invariant"])
        self.assertTrue(rows[(1,)]["d5_singleton_invariant"])
        self.assertTrue(rows[(1, 11)]["highest_label_rule_is_well_defined_in_reduced_symmetry"])

    def test_minimal_breaking_certificate(self):
        cert = self.payload["minimal_breaking_certificate"]
        self.assertEqual(cert["full_aut_group"], [1, 5, 7, 11])
        self.assertFalse(cert["full_aut_d5_singleton_invariant"])
        self.assertEqual(cert["sufficient_subgroups"], [[1], [1, 11]])
        self.assertEqual(cert["nontrivial_sufficient_subgroups"], [[1, 11]])
        self.assertEqual(cert["minimal_nontrivial_sufficient_subgroup"], [1, 11])
        self.assertEqual(cert["minimal_nontrivial_sufficient_subgroup_order"], 2)
        self.assertEqual(cert["excluded_swap_units"], [5, 7])
        self.assertIn("excluding the swap units 5 and 7", cert["interpretation"])

    def test_selector_stack_status_and_guardrails(self):
        payload = self.payload
        status = payload["selector_stack_status"]
        self.assertFalse(status["works_under_full_aut"])
        self.assertTrue(status["works_under_minimal_nontrivial_subgroup"])
        self.assertEqual(status["selected_class_under_reduced_symmetry"], "fifth_step_d5_step_5_7")
        self.assertIn("not derived", status["still_requires_extra_source"])

        proof = payload["proof_reading"]
        self.assertIn("minimal nontrivial sufficient subgroup", proof["finite_gain"])
        self.assertIn("swap d5 with contiguous", proof["negative_result"])
        self.assertIn("reduces the relevant symmetry to {1,11}", proof["conditional_positive_result"])
        self.assertIn("does not derive", proof["remaining_gap"])

        hard_limits = "\n".join(payload["hard_limits"])
        self.assertIn("No identity K_legacy_ont == K_strict_gate", hard_limits)
        self.assertIn("No theorem derives the d5 teacher/self-record trace", hard_limits)
        self.assertIn("No theorem derives the subgroup {1,11}", hard_limits)
        self.assertIn("not strict-core selector closure", hard_limits)
        self.assertIn("No QW-2191 discharge", hard_limits)
        self.assertIn("No ToE closure", hard_limits)
        self.assertTrue(REPORT_MD.exists())


if __name__ == "__main__":
    unittest.main()
