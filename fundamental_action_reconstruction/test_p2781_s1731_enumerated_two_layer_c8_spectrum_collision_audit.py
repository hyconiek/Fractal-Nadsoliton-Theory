import json
import subprocess
import unittest
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]
ROOT = REPO / "fundamental_action_reconstruction"
SCRIPT = ROOT / "p2781_s1731_enumerated_two_layer_c8_spectrum_collision_audit.py"
OUT = ROOT / "generated" / "p2781_s1731_enumerated_two_layer_c8_spectrum_collision_audit.json"
MD = ROOT / "generated" / "p2781_s1731_enumerated_two_layer_c8_spectrum_collision_audit.md"


class P2781EnumeratedTwoLayerC8SpectrumCollisionAuditTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        subprocess.run(["python", str(SCRIPT)], cwd=ROOT, check=True)
        cls.payload = json.loads(OUT.read_text(encoding="utf-8"))

    def test_status_and_input(self):
        self.assertEqual(self.payload["status"], "P2781_ENUMERATED_TWO_LAYER_C8_SPECTRUM_COLLISION_AUDIT_NO_CLOSURE")
        self.assertEqual(self.payload["input_statuses"]["P2780"], "P2780_NONCIRCULANT_16NODE_FULL_SPECTRUM_STRESS_AUDIT_NO_CLOSURE")
        self.assertIn("two-C8-layer", self.payload["audited_question"])

    def test_exhaustive_local_family_has_no_full_spectrum_collision(self):
        witness = self.payload["two_layer_c8_collision_witness"]
        self.assertEqual(witness["base_labeled_candidate_count"], 19)
        self.assertEqual(witness["enumerated_two_layer_candidate_count"], 28)
        self.assertEqual(witness["expanded_labeled_candidate_count"], 47)
        self.assertEqual(witness["expanded_isomorphism_class_count"], 7)
        self.assertEqual(witness["distinct_laplacian_spectrum_count_after_quotient"], 7)
        self.assertEqual(witness["nonisomorphic_cospectral_collision_count"], 0)
        self.assertTrue(witness["full_spectrum_injective_on_expanded_quotient"])
        self.assertEqual(witness["two_layer_member_count_by_representative"]["two_c8_layers_cross_pm0_pm4"], 4)

    def test_acceptance_is_bounded_and_blocks_closure(self):
        acceptance = self.payload["acceptance_matrix"]
        self.assertTrue(acceptance["facts"]["exhaustive_two_layer_c8_shift_pair_family_enumerated"])
        self.assertTrue(acceptance["facts"]["isomorphism_quotient_performed"])
        self.assertTrue(acceptance["accepted_as_exhaustive_local_family_full_spectrum_uniqueness_witness"])
        self.assertFalse(acceptance["accepted_as_strict_spectral_source_law"])
        self.assertFalse(acceptance["accepted_as_canonical_nadsoliton_geometry_source"])
        self.assertFalse(acceptance["accepted_as_ltotal_or_toe_promotion"])
        self.assertIn("strict_nadsoliton_spectral_source_law_exported", acceptance["missing_criteria"])
        self.assertTrue(all(value is False for value in self.payload["decision"]["negative_export_flags"].values()))

    def test_documentation_and_recommendation(self):
        self.assertIn("canonical enumerator for connected 16-node 4-regular graphs", self.payload["decision"]["next_honest_step"])
        self.assertIn("P2697-P2781", self.payload["decision"]["next_honest_step"])
        self.assertIn("P2781/S1731", MD.read_text(encoding="utf-8"))
        self.assertIn("P2781/S1731", (ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md").read_text(encoding="utf-8"))
        self.assertIn("P2781/S1731", (ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md").read_text(encoding="utf-8"))
        self.assertIn("P2781", (REPO / "AGENTS.md").read_text(encoding="utf-8"))


if __name__ == "__main__":
    unittest.main()
