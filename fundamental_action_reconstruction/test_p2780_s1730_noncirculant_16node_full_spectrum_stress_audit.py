import json
import subprocess
import unittest
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]
ROOT = REPO / "fundamental_action_reconstruction"
SCRIPT = ROOT / "p2780_s1730_noncirculant_16node_full_spectrum_stress_audit.py"
OUT = ROOT / "generated" / "p2780_s1730_noncirculant_16node_full_spectrum_stress_audit.json"
MD = ROOT / "generated" / "p2780_s1730_noncirculant_16node_full_spectrum_stress_audit.md"


class P2780Noncirculant16NodeFullSpectrumStressAuditTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        subprocess.run(["python", str(SCRIPT)], cwd=ROOT, check=True)
        cls.payload = json.loads(OUT.read_text(encoding="utf-8"))

    def test_status_and_input(self):
        self.assertEqual(self.payload["status"], "P2780_NONCIRCULANT_16NODE_FULL_SPECTRUM_STRESS_AUDIT_NO_CLOSURE")
        self.assertEqual(self.payload["input_statuses"]["P2779"], "P2779_16NODE_CIRCULANT_FULL_SPECTRUM_QUOTIENT_AUDIT_NO_CLOSURE")
        self.assertIn("non-circulant", self.payload["audited_question"])

    def test_expanded_class_remains_full_spectrum_injective(self):
        witness = self.payload["noncirculant_full_spectrum_stress_witness"]
        self.assertEqual(witness["base_labeled_candidate_count"], 19)
        self.assertEqual(witness["added_labeled_candidate_count"], 4)
        self.assertEqual(witness["expanded_labeled_candidate_count"], 23)
        self.assertEqual(witness["expanded_isomorphism_class_count"], 7)
        self.assertEqual(witness["distinct_laplacian_spectrum_count_after_quotient"], 7)
        self.assertEqual(witness["nonisomorphic_cospectral_collision_count"], 0)
        self.assertTrue(witness["full_spectrum_injective_on_expanded_quotient"])
        self.assertTrue(witness["added_family_spectrum_new_relative_to_p2779_base"])

    def test_acceptance_is_bounded_and_blocks_closure(self):
        acceptance = self.payload["acceptance_matrix"]
        self.assertTrue(acceptance["facts"]["one_broader_16node_4regular_family_added"])
        self.assertTrue(acceptance["facts"]["isomorphism_quotient_performed"])
        self.assertTrue(acceptance["accepted_as_expanded_declared_class_full_spectrum_uniqueness_witness"])
        self.assertFalse(acceptance["accepted_as_strict_spectral_source_law"])
        self.assertFalse(acceptance["accepted_as_canonical_nadsoliton_geometry_source"])
        self.assertFalse(acceptance["accepted_as_ltotal_or_toe_promotion"])
        self.assertIn("strict_nadsoliton_spectral_source_law_exported", acceptance["missing_criteria"])
        self.assertTrue(all(value is False for value in self.payload["decision"]["negative_export_flags"].values()))

    def test_documentation_and_recommendation(self):
        self.assertIn("broader enumerated 16-node 4-regular", self.payload["decision"]["next_honest_step"])
        self.assertIn("P2697-P2780", self.payload["decision"]["next_honest_step"])
        self.assertIn("P2780/S1730", MD.read_text(encoding="utf-8"))
        self.assertIn("P2780/S1730", (ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md").read_text(encoding="utf-8"))
        self.assertIn("P2780/S1730", (ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md").read_text(encoding="utf-8"))
        self.assertIn("P2780", (REPO / "AGENTS.md").read_text(encoding="utf-8"))


if __name__ == "__main__":
    unittest.main()
