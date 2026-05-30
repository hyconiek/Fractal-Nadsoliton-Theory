import json
import subprocess
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "fundamental_action_reconstruction" / "scratch" / "bridge_strict_alpha_nonhistogram_dihedral_quotient_chi11_audit_probe.py"
REPORT_JSON = ROOT / "fundamental_action_reconstruction" / "scratch" / "bridge_strict_alpha_nonhistogram_dihedral_quotient_chi11_audit_report.json"
REPORT_MD = ROOT / "fundamental_action_reconstruction" / "scratch" / "bridge_strict_alpha_nonhistogram_dihedral_quotient_chi11_audit_report.md"


class StrictAlphaNonhistogramDihedralQuotientChi11AuditProbeTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        subprocess.run([sys.executable, str(SCRIPT)], check=True, cwd=ROOT)
        cls.payload = json.loads(REPORT_JSON.read_text(encoding="utf-8"))

    def test_result_kind_and_finite_model(self):
        payload = self.payload
        self.assertEqual(
            payload["result_kind"],
            "SCRATCH_STRICT_ALPHA_NONHISTOGRAM_DIHEDRAL_QUOTIENT_CHI11_AUDIT__BOUNDARY_NOT_A_THEOREM",
        )
        self.assertEqual(
            payload["status"],
            "dihedral-nonhistogram-features-can-carry-chi11-only-with-reduced-symmetry-premise",
        )
        model = payload["finite_model"]
        self.assertEqual(model["ring"], "Z_12")
        self.assertEqual(model["active_count"], 5)
        self.assertEqual(model["support_count"], 792)
        self.assertEqual(model["automorphism_units"], [1, 5, 7, 11])
        self.assertEqual(model["dihedral_units"], [1, 11])
        self.assertEqual(model["dihedral_orbit_count"], 38)
        self.assertEqual(model["full_affine_aut_orbit_count"], 25)
        self.assertEqual(model["dihedral_orbit_size_counts"], {"12": 10, "24": 28})
        self.assertEqual(model["full_affine_orbit_size_counts"], {"12": 4, "24": 11, "48": 10})
        self.assertEqual(model["target_q_power"], "256/243")
        self.assertEqual(model["target_eta"], "9/5")

    def test_branch_orbit_rows_and_certificate(self):
        rows = {row["name"]: row for row in self.payload["branch_orbit_rows"]}
        self.assertEqual(set(rows), {"A1_k1", "A5_k5", "A7_k7", "A11_k11"})
        self.assertEqual(rows["A1_k1"]["support"], [0, 1, 2, 3, 4])
        self.assertEqual(rows["A5_k5"]["support"], [0, 3, 5, 8, 10])
        self.assertEqual(rows["A7_k7"]["support"], [0, 2, 4, 7, 9])
        self.assertEqual(rows["A11_k11"]["support"], [0, 8, 9, 10, 11])
        self.assertEqual(rows["A1_k1"]["required_chi11_value"], -1)
        self.assertEqual(rows["A11_k11"]["required_chi11_value"], -1)
        self.assertEqual(rows["A5_k5"]["required_chi11_value"], 1)
        self.assertEqual(rows["A7_k7"]["required_chi11_value"], 1)
        self.assertEqual(rows["A1_k1"]["dihedral_orbit_index"], rows["A11_k11"]["dihedral_orbit_index"])
        self.assertEqual(rows["A5_k5"]["dihedral_orbit_index"], rows["A7_k7"]["dihedral_orbit_index"])
        self.assertNotEqual(rows["A1_k1"]["dihedral_orbit_index"], rows["A5_k5"]["dihedral_orbit_index"])
        self.assertEqual({row["full_affine_orbit_index"] for row in rows.values()}, {0})

        cert = self.payload["branch_full_orbit_certificate"]
        self.assertEqual(cert["full_affine_orbit_index"], 0)
        self.assertEqual(cert["dihedral_component_indices"], [0, 37])
        self.assertEqual(cert["component_gap_necklaces"], [[1, 1, 1, 1, 8], [2, 2, 3, 2, 3]])
        self.assertTrue(cert["unit5_maps_A1_dihedral_component_to_A5_component"])
        self.assertTrue(cert["unit7_maps_A1_dihedral_component_to_A7_component"])
        self.assertTrue(cert["unit11_maps_A1_to_contiguous_reversal_component"])
        self.assertEqual(cert["full_aut_invariant_singleton_A5_not_A1_classifier_count"], 0)

    def test_quotient_summary_and_action_rows(self):
        summary = self.payload["quotient_summary"]
        self.assertEqual(summary["dihedral_boolean_classifier_count_power"], 38)
        self.assertEqual(summary["dihedral_boolean_classifier_total_exact"], "274877906944")
        self.assertEqual(summary["dihedral_branch_chi11_classifier_count_power"], 36)
        self.assertEqual(summary["dihedral_branch_chi11_classifier_total_exact"], "68719476736")
        self.assertEqual(summary["full_aut_boolean_classifier_count_power"], 25)
        self.assertEqual(summary["full_aut_boolean_classifier_total_exact"], "33554432")
        self.assertEqual(summary["full_aut_singleton_A5_not_A1_classifier_count"], 0)
        self.assertEqual(summary["full_aut_blocks_with_two_dihedral_components"], 13)
        self.assertEqual(summary["full_aut_blocks_fixed_on_dihedral_quotient"], 12)
        self.assertEqual(summary["global_integer_chi11_antiinvariant_dimension_on_dihedral_quotient"], 13)
        self.assertEqual(summary["global_nonzero_pm1_chi11_character_count"], 0)

        action_rows = self.payload["dihedral_orbit_action_rows"]
        self.assertEqual(len(action_rows), 38)
        self.assertEqual(action_rows[0]["unit_image_dihedral_orbit_index"], {"1": 0, "5": 37, "7": 37, "11": 0})
        self.assertFalse(action_rows[0]["fixed_by_unit5_mod_dihedral"])
        fixed_count = sum(1 for row in action_rows if row["fixed_by_unit5_mod_dihedral"])
        self.assertEqual(fixed_count, 12)

        decomposition_rows = self.payload["full_affine_orbit_decomposition_rows"]
        self.assertEqual(len(decomposition_rows), 25)
        branch_rows = [row for row in decomposition_rows if row["contains_branch_orbit"]]
        self.assertEqual(len(branch_rows), 1)
        self.assertEqual(branch_rows[0]["dihedral_component_count"], 2)

    def test_witness_interpretation_and_guardrails(self):
        witness = self.payload["branch_compatible_dihedral_sign_witness"]
        self.assertEqual(witness["values_by_branch"], {"A1_k1": -1, "A5_k5": 1, "A7_k7": 1, "A11_k11": -1})
        self.assertIn("unit-axis bit", witness["requires_imported_premise"])

        proof = self.payload["exact_proof_certificate"]
        self.assertIn("C(12,5)=792", proof["finite_domain"])
        self.assertIn("D12 quotient has 38 orbits", proof["quotient_split"])
        self.assertIn("unit 5 swaps them", proof["branch_block"])
        self.assertIn("F(A1)=F(A5)", proof["full_aut_no_go"])
        self.assertIn("matching chi_11", proof["dihedral_positive_boundary"])
        self.assertIn("unit-axis premise", proof["import_diagnosis"])
        self.assertIn("dimension 13", proof["global_character_space"])
        self.assertIn("not all possible strict nadsoliton", proof["scope_limit"])

        interpretation = self.payload["interpretation"]
        self.assertIn("branch-level chi_11 carrier", interpretation["honest_positive"])
        self.assertIn("not full-Aut provenance", interpretation["honest_negative"])
        self.assertIn("beyond histogram-only data", interpretation["relation_to_previous_probe"])

        ontology = self.payload["ontology_guardrail"]
        self.assertIn("nadsoliton itself", ontology["allowed_reading"])
        self.assertIn("No separate informational layer", ontology["forbidden_reading"])

        hard_limits = "\n".join(self.payload["hard_limits"])
        self.assertIn("No identity K_legacy_ont == K_strict_gate", hard_limits)
        self.assertIn("No theorem derives chi_11", hard_limits)
        self.assertIn("No QW-2191 discharge", hard_limits)
        self.assertIn("No ToE closure", hard_limits)
        self.assertTrue(REPORT_MD.exists())


if __name__ == "__main__":
    unittest.main()
