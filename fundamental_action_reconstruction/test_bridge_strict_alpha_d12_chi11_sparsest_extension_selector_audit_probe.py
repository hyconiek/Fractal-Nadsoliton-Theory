import json
import subprocess
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "fundamental_action_reconstruction" / "scratch" / "bridge_strict_alpha_d12_chi11_sparsest_extension_selector_audit_probe.py"
REPORT_JSON = ROOT / "fundamental_action_reconstruction" / "scratch" / "bridge_strict_alpha_d12_chi11_sparsest_extension_selector_audit_report.json"
REPORT_MD = ROOT / "fundamental_action_reconstruction" / "scratch" / "bridge_strict_alpha_d12_chi11_sparsest_extension_selector_audit_report.md"


class StrictAlphaD12Chi11SparsestExtensionSelectorAuditProbeTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        subprocess.run([sys.executable, str(SCRIPT)], check=True, cwd=ROOT)
        cls.payload = json.loads(REPORT_JSON.read_text(encoding="utf-8"))

    def test_result_kind_and_finite_model(self):
        payload = self.payload
        self.assertEqual(
            payload["result_kind"],
            "SCRATCH_STRICT_ALPHA_D12_CHI11_SPARSEST_EXTENSION_SELECTOR_AUDIT__CONDITIONAL_NOT_A_THEOREM",
        )
        self.assertEqual(
            payload["status"],
            "sparsest-branch-normalized-d12-chi11-extension-is-unique-but-sparsity-is-extra",
        )
        model = payload["finite_model"]
        self.assertEqual(model["ring"], "Z_12")
        self.assertEqual(model["active_count"], 5)
        self.assertEqual(model["support_count"], 792)
        self.assertEqual(model["dihedral_units"], [1, 11])
        self.assertEqual(model["automorphism_units"], [1, 5, 7, 11])
        self.assertEqual(model["d12_orbit_count"], 38)
        self.assertEqual(model["chi11_two_cycle_count"], 13)
        self.assertEqual(model["target_q_power"], "256/243")
        self.assertEqual(model["target_eta"], "9/5")

    def test_sparsity_summary_and_distribution(self):
        summary = self.payload["sparsity_summary"]
        self.assertEqual(summary["branch_cycle_orbit_pair"], [0, 37])
        self.assertEqual(summary["branch_cycle_coefficient_fixed_to"], 1)
        self.assertEqual(summary["free_nonbranch_cycle_count"], 12)
        self.assertEqual(summary["branch_normalized_ternary_function_count"], 531441)
        self.assertEqual(summary["minimum_d12_component_support_size"], 2)
        self.assertEqual(summary["minimum_support_function_count"], 1)
        self.assertTrue(summary["unique_minimum_is_branch_generator"])
        self.assertEqual(summary["maximum_d12_component_support_size"], 26)

        distribution = self.payload["support_size_distribution"]
        self.assertEqual(len(distribution), 13)
        self.assertEqual(sum(row["function_count"] for row in distribution), 531441)
        self.assertEqual(distribution[0], {"nonzero_free_coefficients": 0, "d12_component_support_size": 2, "function_count": 1})
        self.assertEqual(distribution[1], {"nonzero_free_coefficients": 1, "d12_component_support_size": 4, "function_count": 24})
        self.assertEqual(distribution[2], {"nonzero_free_coefficients": 2, "d12_component_support_size": 6, "function_count": 264})
        self.assertEqual(distribution[-1], {"nonzero_free_coefficients": 12, "d12_component_support_size": 26, "function_count": 4096})

    def test_minimum_witness_and_branch_rows(self):
        witness = self.payload["minimum_support_witness"]
        self.assertEqual(witness["nonzero_cycle_coefficients"], {"branch_cycle_[0,37]": 1})
        self.assertEqual(witness["values_by_d12_orbit"], {"0": -1, "37": 1})
        self.assertEqual(witness["values_by_branch"], {"A1_k1": -1, "A5_k5": 1, "A7_k7": 1, "A11_k11": -1})
        self.assertIn("sparsest-extension selector", witness["requires_imported_premise"])

        branches = {row["name"]: row for row in self.payload["branch_rows"]}
        self.assertEqual(branches["A1_k1"]["dihedral_orbit_index"], 0)
        self.assertEqual(branches["A11_k11"]["dihedral_orbit_index"], 0)
        self.assertEqual(branches["A5_k5"]["dihedral_orbit_index"], 37)
        self.assertEqual(branches["A7_k7"]["dihedral_orbit_index"], 37)
        self.assertEqual(branches["A1_k1"]["required_chi11_value"], -1)
        self.assertEqual(branches["A5_k5"]["required_chi11_value"], 1)

    def test_proof_interpretation_and_guardrails(self):
        proof = self.payload["exact_proof_certificate"]
        self.assertIn("C(12,5)=792", proof["finite_domain"])
        self.assertIn("13 two-cycle basis", proof["module_reuse_boundary"])
        self.assertIn("coefficient to +1", proof["branch_constraint"])
        self.assertIn("3^12=531441", proof["free_coefficients"])
        self.assertIn("2+2k", proof["support_formula"])
        self.assertIn("unique minimum-support extension", proof["unique_minimum"])
        self.assertIn("sparsest-extension rule is an added selector", proof["strict_limit"])

        interpretation = self.payload["interpretation"]
        self.assertIn("uniquely selected", interpretation["honest_positive"])
        self.assertIn("not derived from strict geometry", interpretation["honest_negative"])
        self.assertIn("support-size distribution", interpretation["relation_to_previous_probe"])

        ontology = self.payload["ontology_guardrail"]
        self.assertIn("nadsoliton itself", ontology["allowed_reading"])
        self.assertIn("No separate informational layer", ontology["forbidden_reading"])

        hard_limits = "\n".join(self.payload["hard_limits"])
        self.assertIn("No identity K_legacy_ont == K_strict_gate", hard_limits)
        self.assertIn("No theorem derives chi_11", hard_limits)
        self.assertIn("No theorem derives a sparsest-extension selector", hard_limits)
        self.assertIn("No QW-2191 discharge", hard_limits)
        self.assertIn("No ToE closure", hard_limits)
        self.assertTrue(REPORT_MD.exists())


if __name__ == "__main__":
    unittest.main()
