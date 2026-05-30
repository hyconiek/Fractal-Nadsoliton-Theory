from __future__ import annotations

import json
import subprocess
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "fundamental_action_reconstruction" / "scratch" / "bridge_strict_alpha_reynolds_annihilator_chi11_matrix_certificate_probe.py"
REPORT_JSON = ROOT / "fundamental_action_reconstruction" / "scratch" / "bridge_strict_alpha_reynolds_annihilator_chi11_matrix_certificate_report.json"
REPORT_MD = ROOT / "fundamental_action_reconstruction" / "scratch" / "bridge_strict_alpha_reynolds_annihilator_chi11_matrix_certificate_report.md"


class ReynoldsAnnihilatorChi11MatrixCertificateTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        subprocess.run([sys.executable, str(SCRIPT)], check=True, cwd=ROOT)
        cls.payload = json.loads(REPORT_JSON.read_text(encoding="utf-8"))

    def test_result_kind_and_finite_model(self):
        payload = self.payload
        self.assertEqual(
            payload["result_kind"],
            "SCRATCH_STRICT_ALPHA_REYNOLDS_ANNIHILATOR_CHI11_MATRIX_CERTIFICATE__NO_FALSE_PASS",
        )
        self.assertEqual(payload["status"], "exact-integer-matrix-annihilator-computed-on-translation-quotient")
        model = payload["finite_model"]
        self.assertEqual(model["ring"], "Z_12")
        self.assertEqual(model["active_count"], 5)
        self.assertEqual(model["support_count"], 792)
        self.assertEqual(model["translation_orbit_count"], 66)
        self.assertEqual(model["matrix_size"], 66)
        self.assertEqual(model["residual_unit_group"], [1, 5, 7, 11])
        self.assertEqual(model["residual_unit_orbit_count"], 25)
        self.assertEqual(model["residual_unit_orbit_size_counts"], {"1": 4, "2": 11, "4": 10})
        self.assertEqual(model["projector_denominator"], 4)
        self.assertEqual(model["target_q_power"], "256/243")
        self.assertEqual(model["target_eta"], "9/5")

    def test_matrix_annihilator_certificate(self):
        matrix = self.payload["matrix_certificate"]
        self.assertEqual(matrix["reynolds_numerator_rank"], 25)
        self.assertEqual(matrix["chi11_numerator_rank"], 13)
        self.assertTrue(matrix["reynolds_times_chi11_numerator_is_zero_matrix"])
        self.assertTrue(matrix["chi11_times_reynolds_numerator_is_zero_matrix"])
        self.assertEqual(matrix["annihilator_zero_entry_count"], 66 * 66)
        self.assertEqual(matrix["annihilator_total_entry_count"], 66 * 66)
        self.assertIn("before dividing by 4", matrix["identity"])

    def test_branch_generator_certificate(self):
        branch = self.payload["branch_generator_certificate"]
        self.assertEqual(branch["A1_translation_orbit_index"], 0)
        self.assertEqual(branch["A11_translation_orbit_index"], 0)
        self.assertEqual(branch["A5_translation_orbit_index"], 65)
        self.assertEqual(branch["A7_translation_orbit_index"], 65)
        self.assertEqual(branch["quotient_pair_indices"], [0, 65])
        self.assertEqual(branch["normalized_values_by_branch"], {"A1_k1": -1, "A11_k11": -1, "A5_k5": 1, "A7_k7": 1})
        self.assertEqual(branch["compact_translation_vector_nonzero_entries"], {"0": -1, "65": 1})
        self.assertEqual(branch["A1_representative_support"], [0, 1, 2, 3, 4])
        self.assertEqual(branch["A5_representative_support"], [0, 2, 4, 7, 9])
        self.assertTrue(branch["reynolds_numerator_on_branch_is_zero"])
        self.assertTrue(branch["chi11_numerator_on_branch_equals_four_times_branch"])
        self.assertEqual(branch["reynolds_numerator_image_nonzero_entries"], {})
        self.assertEqual(branch["chi11_numerator_image_nonzero_entries"], {"0": -4, "65": 4})
        self.assertTrue(all(row["branch_sum"] == 0 for row in branch["full_unit_orbit_branch_sums"]))

    def test_chi11_basis_rows_proof_guardrails_and_markdown(self):
        basis = self.payload["chi11_basis_annihilation_summary"]
        self.assertEqual(basis["chi11_basis_row_count"], 13)
        self.assertTrue(basis["all_chi11_basis_rows_have_zero_reynolds_sum"])
        self.assertTrue(all(row["unit_orbit_signed_sum"] == 0 for row in basis["basis_rows"]))
        self.assertEqual(basis["basis_rows"][0]["translation_orbit_indices"], [0, 65])
        self.assertEqual(basis["basis_rows"][0]["values_by_translation_orbit_index"], {"0": 1, "65": -1})

        proof = self.payload["exact_proof_certificate"]
        self.assertIn("C(12,5)=792", proof["finite_domain"])
        self.assertIn("66x66 permutation matrix", proof["matrix_construction"])
        self.assertIn("zero matrices", proof["annihilator_identity"])
        self.assertIn("R_num*v=0", proof["branch_witness"])
        self.assertIn("extra premise", proof["logical_boundary"])

        interpretation = self.payload["interpretation"]
        self.assertIn("explicit integer matrix annihilator", interpretation["honest_positive"])
        self.assertIn("loses the branch polarity", interpretation["honest_negative"])
        self.assertIn("does not repeat subgroup", interpretation["relation_to_previous_probe"])

        ontology = self.payload["ontology_guardrail"]
        self.assertIn("nadsoliton itself", ontology["allowed_reading"])
        self.assertIn("No separate informational layer", ontology["forbidden_reading"])

        hard_limits = "\n".join(self.payload["hard_limits"])
        self.assertIn("No identity K_legacy_ont == K_strict_gate", hard_limits)
        self.assertIn("No theorem derives chi_11", hard_limits)
        self.assertIn("No theorem derives a full-Aut internal chi_11 polarity source", hard_limits)
        self.assertIn("No QW-2191 discharge", hard_limits)
        self.assertIn("No ToE closure", hard_limits)
        self.assertTrue(REPORT_MD.exists())


if __name__ == "__main__":
    unittest.main()
