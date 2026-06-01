import json
import subprocess
import sys
import unittest
from pathlib import Path


class StrictCompletionToEProperSubsetObstructionCertificateProbeTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.repo = Path(__file__).resolve().parents[1]
        cls.script = cls.repo / "fundamental_action_reconstruction" / "scratch" / "bridge_strict_completion_toe_proper_subset_obstruction_certificate_probe.py"
        subprocess.run([sys.executable, str(cls.script)], check=True, cwd=cls.repo)
        cls.report = cls.repo / "fundamental_action_reconstruction" / "scratch" / "bridge_strict_completion_toe_proper_subset_obstruction_certificate_report.json"
        cls.markdown = cls.repo / "fundamental_action_reconstruction" / "scratch" / "bridge_strict_completion_toe_proper_subset_obstruction_certificate_report.md"
        cls.payload = json.loads(cls.report.read_text(encoding="utf-8"))

    def test_result_kind_sources_and_counts(self):
        payload = self.payload
        self.assertEqual(
            payload["result_kind"],
            "SCRATCH_STRICT_COMPLETION_TOE_PROPER_SUBSET_OBSTRUCTION_CERTIFICATE__ALL_127_PROPER_SUBSETS_FAIL_NO_CLOSURE",
        )
        self.assertTrue(self.report.exists())
        self.assertTrue(self.markdown.exists())
        self.assertEqual(
            set(payload["source_reports"]),
            {"frontier_truth_table", "frontier_cut", "target_signature_lattice", "toe_potential_readiness"},
        )
        self.assertEqual(payload["proper_subset_count"], 127)
        self.assertEqual(payload["proper_subset_count_by_weight"], {"0": 1, "1": 7, "2": 21, "3": 35, "4": 35, "5": 21, "6": 7})
        self.assertEqual(payload["max_proper_true_targets"], 2)
        self.assertEqual(payload["max_proper_target_signatures"], ["0110", "1010"])
        self.assertEqual(len(payload["nearest_miss_rows"]), 7)

    def test_cross_checks_and_summary(self):
        payload = self.payload
        self.assertEqual(payload["status"], "PASS")
        self.assertTrue(payload["all_cross_checks_pass"])
        self.assertTrue(all(payload["cross_checks"].values()))
        checks = payload["cross_checks"]
        self.assertTrue(checks["source_reports_present"])
        self.assertTrue(checks["truth_table_count_pass"])
        self.assertTrue(checks["proper_subset_count_pass"])
        self.assertTrue(checks["full_set_only_toe_pass"])
        self.assertTrue(checks["nearest_miss_count_pass"])
        self.assertTrue(checks["nearest_miss_rows_fail_toe"])
        self.assertTrue(checks["missing_atom_rows_cover_frontier"])
        self.assertTrue(checks["minimal_cut_matches_truth_table"])
        self.assertTrue(checks["target_lattice_full_signature_only"])
        self.assertTrue(checks["readiness_certificate_inherited"])
        self.assertTrue(checks["toe_audit_doc_mentions_joint_requirement"])
        self.assertTrue(checks["hard_limits_preserved"])

        summary = payload["proper_subset_obstruction_summary"]
        self.assertTrue(summary["all_127_proper_subsets_fail_toe"])
        self.assertEqual(summary["nearest_miss_count"], 7)
        self.assertTrue(summary["six_atom_packages_fail_toe"])
        self.assertTrue(summary["all_seven_atoms_required"])
        self.assertEqual(summary["max_proper_true_targets"], 2)
        self.assertFalse(summary["toe_closure_claimed"])

    def test_nearest_miss_explanations_and_limits(self):
        payload = self.payload
        rows = {row["missing_atom"]: row for row in payload["nearest_miss_rows"]}
        self.assertEqual(rows["chi11_selector_source"]["target_signature"], "1000")
        self.assertFalse(rows["chi11_selector_source"]["selector_qw2191_closure"])
        self.assertIn("selector and role atom missing", rows["chi11_selector_source"]["deficit_explanation"])
        self.assertEqual(rows["strict_dynamical_source_for_A_P_D"]["target_signature"], "0110")
        self.assertFalse(rows["strict_dynamical_source_for_A_P_D"]["bridge_theorem_level_closure"])
        self.assertEqual(rows["alpha_geo_electroweak_role_theorem"]["target_signature"], "1010")
        self.assertFalse(rows["alpha_geo_electroweak_role_theorem"]["role_transfer_theorem_level_closure"])

        proof = payload["proof_certificate"]
        self.assertIn("all 127 proper subsets", proof["enumeration_step"])
        self.assertIn("seven six-atom nearest misses", proof["nearest_miss_step"])
        self.assertIn("No proper subset reaches ToE signature 1111", proof["signature_step"])
        self.assertIn("exports no new theorem atom", proof["limit_step"])
        hard_limits = "\n".join(payload["hard_limits"])
        self.assertIn("No proper subset is promoted to ToE closure", hard_limits)
        self.assertIn("No new theorem atom", hard_limits)
        self.assertIn("No QW-2191 selector discharge", hard_limits)
        self.assertIn("No ToE closure", hard_limits)


if __name__ == "__main__":
    unittest.main()
