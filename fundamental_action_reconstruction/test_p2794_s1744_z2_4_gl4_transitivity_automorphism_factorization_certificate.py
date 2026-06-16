import json
import unittest
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]
ROOT = REPO / "fundamental_action_reconstruction"
GEN = ROOT / "generated"
JSON_PATH = GEN / "p2794_s1744_z2_4_gl4_transitivity_automorphism_factorization_certificate.json"
MD_PATH = GEN / "p2794_s1744_z2_4_gl4_transitivity_automorphism_factorization_certificate.md"


class P2794Z2GL4TransitivityAutomorphismFactorizationCertificateTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        if not JSON_PATH.exists():
            import subprocess
            subprocess.run(["python", str(ROOT / "p2794_s1744_z2_4_gl4_transitivity_automorphism_factorization_certificate.py")], check=True, cwd=ROOT)
        cls.payload = json.loads(JSON_PATH.read_text(encoding="utf-8"))
        cls.witness = cls.payload["gl4_transitivity_witness"]
        cls.acceptance = cls.payload["acceptance_matrix"]

    def test_status_and_input(self):
        self.assertEqual(self.payload["status"], "P2794_Z2_4_GL4_TRANSITIVITY_AUTOMORPHISM_FACTORIZATION_CERTIFICATE_NO_CLOSURE")
        self.assertEqual(self.payload["input_statuses"]["P2793"], "P2793_Z2_4_FOUR_GENERATOR_CAYLEY_SUBCLASS_EXHAUSTION_CERTIFICATE_NO_CLOSURE")

    def test_gl4_order_and_basis_transitivity(self):
        self.assertEqual(self.witness["all_binary_matrix_count"], 65536)
        self.assertEqual(self.witness["gl4_2_order"], 20160)
        self.assertEqual(self.witness["gl4_2_order_formula_factors"], [15, 14, 12, 8])
        self.assertEqual(self.witness["gl4_2_order_formula_product"], 20160)
        self.assertEqual(self.witness["unordered_basis_count"], 840)
        self.assertEqual(self.witness["basis_count_formula_value"], 840)
        self.assertEqual(self.witness["preimage_count_histogram"], {"24": 840})
        self.assertTrue(self.witness["all_unordered_bases_have_24_maps_to_standard_basis"])
        self.assertTrue(self.witness["all_constructed_preimage_maps_are_in_gl4"])
        self.assertTrue(self.witness["all_constructed_preimage_maps_hit_standard_basis"])

    def test_automorphism_factorization(self):
        self.assertEqual(self.witness["standard_graph_edge_count"], 32)
        self.assertEqual(self.witness["linear_generator_permutation_count"], 24)
        self.assertEqual(self.witness["translation_count"], 16)
        self.assertEqual(self.witness["affine_automorphism_count_from_factorization"], 384)
        self.assertTrue(self.witness["all_factored_affine_maps_preserve_standard_edges"])
        self.assertTrue(self.witness["matches_p2793_automorphism_group_size"])

    def test_acceptance_blocks_closure(self):
        self.assertTrue(self.acceptance["accepted_as_gl4_transitivity_and_automorphism_factorization_certificate"])
        self.assertFalse(self.acceptance["accepted_as_full_16node_canonical_generator_certificate"])
        self.assertFalse(self.acceptance["accepted_as_strict_spectral_source_law"])
        self.assertFalse(self.acceptance["accepted_as_canonical_nadsoliton_geometry_source"])
        self.assertFalse(self.acceptance["accepted_as_ltotal_or_toe_promotion"])
        self.assertIn("P2697-P2794", self.payload["decision"]["next_honest_step"])
        self.assertTrue(all(value is False for value in self.payload["decision"]["negative_export_flags"].values()))

    def test_documentation_updated(self):
        self.assertIn("P2794/S1744", MD_PATH.read_text(encoding="utf-8"))
        self.assertIn("P2794/S1744", (ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md").read_text(encoding="utf-8"))
        self.assertIn("P2794/S1744", (ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md").read_text(encoding="utf-8"))
        self.assertIn("P2794", (REPO / "AGENTS.md").read_text(encoding="utf-8"))


if __name__ == "__main__":
    unittest.main()
