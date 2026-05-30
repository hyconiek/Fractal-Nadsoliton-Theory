import json
import subprocess
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "fundamental_action_reconstruction" / "scratch" / "bridge_strict_alpha_chi11_premise_dependency_lattice_audit_probe.py"
REPORT_JSON = ROOT / "fundamental_action_reconstruction" / "scratch" / "bridge_strict_alpha_chi11_premise_dependency_lattice_audit_report.json"
REPORT_MD = ROOT / "fundamental_action_reconstruction" / "scratch" / "bridge_strict_alpha_chi11_premise_dependency_lattice_audit_report.md"


class StrictAlphaChi11PremiseDependencyLatticeAuditProbeTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        subprocess.run([sys.executable, str(SCRIPT)], check=True, cwd=ROOT)
        cls.payload = json.loads(REPORT_JSON.read_text(encoding="utf-8"))

    def test_result_kind_and_finite_model(self):
        payload = self.payload
        self.assertEqual(
            payload["result_kind"],
            "SCRATCH_STRICT_ALPHA_CHI11_PREMISE_DEPENDENCY_LATTICE_AUDIT__NO_STRICT_FULL_AUT_SOURCE",
        )
        self.assertEqual(
            payload["status"],
            "minimal-premise-antichain-computed-for-existing-chi11-support-audits",
        )
        model = payload["finite_model"]
        self.assertEqual(model["ring"], "Z_12")
        self.assertEqual(model["active_count"], 5)
        self.assertEqual(model["support_count"], 792)
        self.assertEqual(model["premise_count"], 8)
        self.assertEqual(model["premise_subset_count"], 256)
        self.assertEqual(model["target_q_power"], "256/243")
        self.assertEqual(model["target_eta"], "9/5")

    def test_cross_report_checks(self):
        checks = self.payload["cross_report_checks"]
        self.assertEqual(checks["input_report_count"], 7)
        self.assertTrue(checks["all_input_reports_loaded"])
        self.assertTrue(checks["support_count_agrees"])
        self.assertEqual(checks["full_aut_orbit_count"], 25)
        self.assertEqual(checks["d12_orbit_count"], 38)
        self.assertEqual(checks["cyclic_translation_orbit_count"], 66)
        self.assertTrue(checks["full_aut_kills_chi11"])
        self.assertTrue(checks["unit5_and_unit7_are_killing_additions"])
        self.assertTrue(checks["full_aut_block_amplitude_locates_unique_block"])
        self.assertTrue(checks["full_aut_block_amplitude_exports_no_polarity"])
        self.assertTrue(checks["nonhistogram_singleton_A5_not_A1_count_zero"])
        self.assertEqual(checks["cyclic_chi11_dimension"], 13)
        self.assertEqual(checks["cyclic_full_aut_trivial_intersection_with_chi11_rank"], 0)
        self.assertEqual(checks["d12_chi11_rank"], 13)
        self.assertTrue(checks["d12_full_aut_intersection_rank_zero"])
        self.assertTrue(checks["sparsest_unique_minimum_is_branch_generator"])
        self.assertTrue(checks["max_shell_imbalance_unique_branch_maximizer"])
        self.assertTrue(checks["max_shell_imbalance_requires_shell_label"])

    def test_minimal_premise_antichains(self):
        antichains = self.payload["minimal_premise_antichains"]
        self.assertEqual(
            antichains["locate_branch_full_aut_block_without_polarity"],
            [["full_aut_unoriented_block_amplitude"]],
        )
        self.assertEqual(
            antichains["host_nonzero_cyclic_chi11_character_space"],
            [["cyclic_translation_quotient", "chi11_unit_character_choice"]],
        )
        self.assertEqual(
            antichains["host_d12_chi11_covariant_module"],
            [["d12_reduced_quotient", "chi11_unit_character_choice"]],
        )
        self.assertEqual(
            antichains["branch_normalized_d12_chi11_family"],
            [["d12_reduced_quotient", "chi11_unit_character_choice", "branch_normalization_constraint"]],
        )
        self.assertEqual(
            antichains["unique_branch_generator_by_sparsest_extension"],
            [[
                "d12_reduced_quotient",
                "chi11_unit_character_choice",
                "branch_normalization_constraint",
                "sparsest_extension_selector",
            ]],
        )
        self.assertEqual(
            antichains["unique_branch_generator_by_max_shell_imbalance"],
            [[
                "d12_reduced_quotient",
                "chi11_unit_character_choice",
                "branch_normalization_constraint",
                "shell_label_d1_d5_axis",
                "max_shell_imbalance_selector",
            ]],
        )
        self.assertEqual(antichains["strict_full_aut_internal_chi11_polarity_source"], [])

    def test_rows_proof_guardrails_and_markdown(self):
        rows = {row["outcome"]: row for row in self.payload["outcome_rows"]}
        self.assertIsNone(rows["strict_full_aut_internal_chi11_polarity_source"]["minimal_premise_count"])
        self.assertFalse(rows["strict_full_aut_internal_chi11_polarity_source"]["realized_by_current_reports"])
        self.assertEqual(rows["unique_branch_generator_by_sparsest_extension"]["minimal_premise_count"], 4)
        self.assertEqual(rows["unique_branch_generator_by_max_shell_imbalance"]["minimal_premise_count"], 5)

        proof = self.payload["exact_proof_certificate"]
        self.assertIn("C(12,5)=792", proof["finite_domain"])
        self.assertIn("P2165", proof["not_duplicate_of_general_qw2191_lattice"])
        self.assertIn("exports_chi11_polarity=false", proof["full_aut_no_polarity"])
        self.assertIn("dimension 13", proof["reduced_quotient_boundary"])
        self.assertIn("sparsest extension", proof["selector_boundary"])
        self.assertIn("No row", proof["strict_limit"])

        interpretation = self.payload["interpretation"]
        self.assertIn("finite premise-dependency antichain", interpretation["honest_positive"])
        self.assertIn("no minimal set", interpretation["honest_negative"])
        self.assertIn("not another duplicate enumeration", interpretation["relation_to_previous_probe"])

        ontology = self.payload["ontology_guardrail"]
        self.assertIn("nadsoliton itself", ontology["allowed_reading"])
        self.assertIn("No separate informational layer", ontology["forbidden_reading"])

        hard_limits = "\n".join(self.payload["hard_limits"])
        self.assertIn("No identity K_legacy_ont == K_strict_gate", hard_limits)
        self.assertIn("No theorem derives chi_11", hard_limits)
        self.assertIn("No theorem derives the shell-label", hard_limits)
        self.assertIn("No QW-2191 discharge", hard_limits)
        self.assertIn("No ToE closure", hard_limits)
        self.assertTrue(REPORT_MD.exists())


if __name__ == "__main__":
    unittest.main()
