import json
import subprocess
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "fundamental_action_reconstruction" / "scratch" / "bridge_strict_alpha_cyclic_quotient_character_projection_chi11_certificate_probe.py"
REPORT_JSON = ROOT / "fundamental_action_reconstruction" / "scratch" / "bridge_strict_alpha_cyclic_quotient_character_projection_chi11_certificate_report.json"
REPORT_MD = ROOT / "fundamental_action_reconstruction" / "scratch" / "bridge_strict_alpha_cyclic_quotient_character_projection_chi11_certificate_report.md"


class StrictAlphaCyclicQuotientCharacterProjectionChi11CertificateProbeTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        subprocess.run([sys.executable, str(SCRIPT)], check=True, cwd=ROOT)
        cls.payload = json.loads(REPORT_JSON.read_text(encoding="utf-8"))

    def test_result_kind_and_finite_model(self):
        payload = self.payload
        self.assertEqual(
            payload["result_kind"],
            "SCRATCH_STRICT_ALPHA_CYCLIC_QUOTIENT_CHARACTER_PROJECTION_CHI11_CERTIFICATE__NO_FALSE_PASS",
        )
        self.assertEqual(
            payload["status"],
            "translation-quotient-character-projection-rank-13-for-chi11-but-not-strict-source",
        )
        model = payload["finite_model"]
        self.assertEqual(model["ring"], "Z_12")
        self.assertEqual(model["active_count"], 5)
        self.assertEqual(model["support_count"], 792)
        self.assertEqual(model["translation_orbit_count"], 66)
        self.assertEqual(model["residual_unit_group"], [1, 5, 7, 11])
        self.assertEqual(model["residual_unit_orbit_count"], 25)
        self.assertEqual(model["residual_unit_orbit_size_counts"], {"1": 4, "2": 11, "4": 10})
        self.assertEqual(model["target_q_power"], "256/243")
        self.assertEqual(model["target_eta"], "9/5")

    def test_trace_summary_and_character_dimensions(self):
        trace = self.payload["trace_summary"]
        self.assertEqual(trace["fixed_translation_orbit_counts_by_unit"], {"1": 66, "5": 10, "7": 14, "11": 10})
        self.assertEqual(trace["trivial_full_aut_invariant_dimension"], 25)
        self.assertEqual(trace["chi11_dimension"], 13)
        self.assertEqual(trace["full_aut_trivial_intersection_with_chi11_rank"], 0)

        dims = trace["character_dimension_rows"]
        self.assertEqual(dims["trivial_kernel_full"]["dimension"], 25)
        self.assertEqual(dims["chi11_kernel_{1,11}"]["dimension"], 13)
        self.assertEqual(dims["chi5_kernel_{1,5}"]["dimension"], 13)
        self.assertEqual(dims["chi7_kernel_{1,7}"]["dimension"], 15)
        self.assertEqual(dims["chi11_kernel_{1,11}"]["projection_trace_numerator"], 52)
        self.assertEqual(dims["chi11_kernel_{1,11}"]["values_by_unit"], {"1": 1, "5": -1, "7": -1, "11": 1})

    def test_branch_basis_and_chi11_basis_rows(self):
        cert = self.payload["branch_chi11_basis_certificate"]
        self.assertEqual(cert["basis_index"], 0)
        self.assertTrue(cert["is_branch_basis_block"])
        self.assertEqual(cert["translation_orbit_indices"], [0, 65])
        self.assertEqual(cert["seed_translation_orbit_index"], 0)
        self.assertEqual(cert["seed_representative_support"], [0, 1, 2, 3, 4])
        self.assertEqual(cert["seed_gap_necklace_cyclic"], [1, 1, 1, 1, 8])
        self.assertEqual(cert["stabilizer_units"], [1, 11])
        self.assertEqual(cert["values_by_translation_orbit_index"], {"0": 1, "65": -1})
        self.assertEqual(cert["normalized_values_by_branch"], {"A1_k1": -1, "A5_k5": 1, "A7_k7": 1, "A11_k11": -1})
        self.assertIn("unit-axis", cert["requires_imported_premise"])

        basis = self.payload["chi11_basis_rows"]
        self.assertEqual(len(basis), 13)
        self.assertEqual(sum(1 for row in basis if row["is_branch_basis_block"]), 1)
        self.assertEqual(
            sorted(tuple(row["stabilizer_units"]) for row in basis).count((1, 11)),
            3,
        )

        branches = {row["name"]: row for row in self.payload["branch_rows"]}
        self.assertEqual(branches["A1_k1"]["translation_orbit_index"], 0)
        self.assertEqual(branches["A11_k11"]["translation_orbit_index"], 0)
        self.assertEqual(branches["A5_k5"]["translation_orbit_index"], 65)
        self.assertEqual(branches["A7_k7"]["translation_orbit_index"], 65)

    def test_proof_interpretation_and_guardrails(self):
        proof = self.payload["exact_proof_certificate"]
        self.assertIn("C(12,5)=792", proof["finite_domain"])
        self.assertIn("dim V_chi", proof["character_projection"])
        self.assertIn("rank is 13", proof["chi11_rank"])
        self.assertIn("basis block 0", proof["branch_basis"])
        self.assertIn("intersection", proof["orthogonality_boundary"])
        self.assertIn("does not derive the unit-axis premise", proof["strict_limit"])

        interpretation = self.payload["interpretation"]
        self.assertIn("13-dimensional chi_11 sector", interpretation["honest_positive"])
        self.assertIn("zero intersection", interpretation["honest_negative"])
        self.assertIn("character-projection", interpretation["relation_to_previous_probe"])

        ontology = self.payload["ontology_guardrail"]
        self.assertIn("nadsoliton itself", ontology["allowed_reading"])
        self.assertIn("No separate informational layer", ontology["forbidden_reading"])

        hard_limits = "\n".join(self.payload["hard_limits"])
        self.assertIn("No identity K_legacy_ont == K_strict_gate", hard_limits)
        self.assertIn("No theorem derives chi_11", hard_limits)
        self.assertIn("No QW-2191 discharge", hard_limits)
        self.assertIn("No ToE closure", hard_limits)
        self.assertIn("does not supply an internal strict source", hard_limits)
        self.assertTrue(REPORT_MD.exists())


if __name__ == "__main__":
    unittest.main()
