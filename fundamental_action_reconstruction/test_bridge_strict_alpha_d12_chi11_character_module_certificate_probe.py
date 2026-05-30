import json
import subprocess
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "fundamental_action_reconstruction" / "scratch" / "bridge_strict_alpha_d12_chi11_character_module_certificate_probe.py"
REPORT_JSON = ROOT / "fundamental_action_reconstruction" / "scratch" / "bridge_strict_alpha_d12_chi11_character_module_certificate_report.json"
REPORT_MD = ROOT / "fundamental_action_reconstruction" / "scratch" / "bridge_strict_alpha_d12_chi11_character_module_certificate_report.md"


class StrictAlphaD12Chi11CharacterModuleCertificateProbeTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        subprocess.run([sys.executable, str(SCRIPT)], check=True, cwd=ROOT)
        cls.payload = json.loads(REPORT_JSON.read_text(encoding="utf-8"))

    def test_result_kind_and_finite_model(self):
        payload = self.payload
        self.assertEqual(
            payload["result_kind"],
            "SCRATCH_STRICT_ALPHA_D12_CHI11_CHARACTER_MODULE_CERTIFICATE__BOUNDARY_NOT_A_THEOREM",
        )
        self.assertEqual(
            payload["status"],
            "exact-d12-quotient-chi11-module-rank-13-but-requires-unit-axis-premise",
        )
        model = payload["finite_model"]
        self.assertEqual(model["ring"], "Z_12")
        self.assertEqual(model["active_count"], 5)
        self.assertEqual(model["support_count"], 792)
        self.assertEqual(model["dihedral_units"], [1, 11])
        self.assertEqual(model["automorphism_units"], [1, 5, 7, 11])
        self.assertEqual(model["d12_orbit_count"], 38)
        self.assertEqual(model["residual_unit5_permutation_size"], 38)
        self.assertEqual(model["target_q_power"], "256/243")
        self.assertEqual(model["target_eta"], "9/5")

    def test_module_summary_and_basis_rows(self):
        summary = self.payload["module_summary"]
        self.assertEqual(summary["tau_two_cycle_count"], 13)
        self.assertEqual(summary["tau_fixed_orbit_count"], 12)
        self.assertEqual(summary["integer_chi11_covariant_module_rank"], 13)
        self.assertEqual(summary["full_aut_invariant_intersection_rank"], 0)
        self.assertEqual(summary["global_nonzero_pm1_character_count"], 0)
        self.assertEqual(summary["ternary_minus_zero_plus_chi11_covariant_count"], 1594323)
        self.assertEqual(summary["branch_normalized_ternary_chi11_covariant_count"], 531441)
        self.assertEqual(summary["fixed_orbit_forced_value"], 0)

        basis_rows = self.payload["two_cycle_basis_rows"]
        fixed_rows = self.payload["fixed_orbit_rows"]
        self.assertEqual(len(basis_rows), 13)
        self.assertEqual(len(fixed_rows), 12)
        self.assertEqual([row["cycle_index"] for row in basis_rows], list(range(13)))
        self.assertTrue(all(row["basis_values_low_high"] == [-1, 1] for row in basis_rows))
        self.assertTrue(all(row["forced_chi11_value"] == 0 for row in fixed_rows))
        self.assertEqual(sum(1 for row in basis_rows if row["is_branch_chi11_generator"]), 1)

    def test_branch_generator_certificate(self):
        branches = {row["name"]: row for row in self.payload["branch_rows"]}
        self.assertEqual(branches["A1_k1"]["support"], [0, 1, 2, 3, 4])
        self.assertEqual(branches["A5_k5"]["support"], [0, 3, 5, 8, 10])
        self.assertEqual(branches["A7_k7"]["support"], [0, 2, 4, 7, 9])
        self.assertEqual(branches["A11_k11"]["support"], [0, 8, 9, 10, 11])
        self.assertEqual(branches["A1_k1"]["dihedral_orbit_index"], 0)
        self.assertEqual(branches["A11_k11"]["dihedral_orbit_index"], 0)
        self.assertEqual(branches["A5_k5"]["dihedral_orbit_index"], 37)
        self.assertEqual(branches["A7_k7"]["dihedral_orbit_index"], 37)

        cert = self.payload["branch_generator_certificate"]
        self.assertEqual(cert["cycle_index"], 0)
        self.assertEqual(cert["orbit_pair"], [0, 37])
        self.assertEqual(cert["basis_values_low_high"], [-1, 1])
        self.assertEqual(cert["low_gap_necklace"], [1, 1, 1, 1, 8])
        self.assertEqual(cert["high_gap_necklace"], [2, 2, 3, 2, 3])
        self.assertEqual(cert["values_by_branch"], {"A1_k1": -1, "A5_k5": 1, "A7_k7": 1, "A11_k11": -1})
        self.assertIn("unit-axis orientation", cert["requires_imported_premise"])
        self.assertIn("not full-Aut", cert["requires_imported_premise"])

    def test_proof_interpretation_and_guardrails(self):
        proof = self.payload["exact_proof_certificate"]
        self.assertIn("C(12,5)=792", proof["finite_domain"])
        self.assertIn("units 5 and 7 induce the same involution", proof["residual_action"])
        self.assertIn("f(tau(O))=-f(O)", proof["module_equations"])
        self.assertIn("rank 13", proof["rank_certificate"])
        self.assertIn("f(tau(O))=f(O)", proof["full_aut_intersection"])
        self.assertIn("A1/A11 orbit and A5/A7 orbit", proof["branch_generator"])
        self.assertIn("unit-axis orientation", proof["not_strict_source"])

        interpretation = self.payload["interpretation"]
        self.assertIn("complete D12-quotient", interpretation["honest_positive"])
        self.assertIn("full-Aut invariant support data", interpretation["honest_negative"])
        self.assertIn("solves the full character module", interpretation["relation_to_previous_probe"])

        ontology = self.payload["ontology_guardrail"]
        self.assertIn("nadsoliton itself", ontology["allowed_reading"])
        self.assertIn("No separate informational layer", ontology["forbidden_reading"])

        hard_limits = "\n".join(self.payload["hard_limits"])
        self.assertIn("No identity K_legacy_ont == K_strict_gate", hard_limits)
        self.assertIn("No theorem derives chi_11", hard_limits)
        self.assertIn("No QW-2191 discharge", hard_limits)
        self.assertIn("No ToE closure", hard_limits)
        self.assertIn("not a complete strict-source provenance theorem", hard_limits)
        self.assertTrue(REPORT_MD.exists())


if __name__ == "__main__":
    unittest.main()
