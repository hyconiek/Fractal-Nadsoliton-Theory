import json
import subprocess
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
JSON_PATH = ROOT / "generated" / "p2805_s1755_girth4_spectral_collision_isomorphism_refinement.json"
MD_PATH = ROOT / "generated" / "p2805_s1755_girth4_spectral_collision_isomorphism_refinement.md"


class P2805Girth4SpectralCollisionIsomorphismRefinementTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        subprocess.run([sys.executable, str(ROOT / "p2805_s1755_girth4_spectral_collision_isomorphism_refinement.py")], check=True, cwd=ROOT)
        cls.payload = json.loads(JSON_PATH.read_text(encoding="utf-8"))
        cls.refinement = cls.payload["spectral_collision_isomorphism_refinement"]
        cls.acceptance = cls.payload["acceptance_matrix"]

    def test_status_and_p2804_input(self):
        self.assertEqual(
            self.payload["status"],
            "P2805_GIRTH4_SPECTRAL_COLLISION_ISOMORPHISM_REFINEMENT_NO_CANONICAL_LABELS_NO_SOURCE_LAW_NO_CLOSURE",
        )
        self.assertEqual(
            self.payload["input_statuses"]["P2804"],
            "P2804_GIRTH4_SPECTRAL_COMPLEMENT_QUOTIENT_AUDIT_NO_ISOMORPHISM_NO_SOURCE_LAW_NO_CLOSURE",
        )
        self.assertTrue(self.refinement["p2804_accepts_spectral_complement_quotient"])

    def test_refines_all_spectral_collision_classes(self):
        self.assertEqual(self.refinement["decoded_graph_count"], 16828)
        self.assertEqual(self.refinement["spectral_pair_class_count"], 16211)
        self.assertEqual(self.refinement["spectral_pair_singleton_class_count"], 15633)
        self.assertEqual(self.refinement["spectral_pair_collision_class_count"], 578)
        self.assertEqual(self.refinement["spectral_pair_collision_graph_count"], 1195)
        self.assertEqual(self.refinement["isomorphism_pairwise_checks_against_component_representatives"], 660)
        self.assertEqual(self.refinement["negative_isomorphism_rejections_inside_collisions"], 660)

    def test_no_duplicate_isomorphic_records_found(self):
        self.assertEqual(self.refinement["positive_isomorphism_matches_inside_collisions"], 0)
        self.assertEqual(self.refinement["duplicate_isomorphism_component_count"], 0)
        self.assertEqual(self.refinement["component_size_histogram_inside_collisions"], {"1": 1195})
        self.assertEqual(self.refinement["resolved_isomorphism_classes_inside_spectral_collisions"], 1195)
        self.assertEqual(self.refinement["resolved_total_isomorphism_classes_after_refinement"], 16828)

    def test_acceptance_boundaries(self):
        self.assertTrue(self.acceptance["accepted_as_spectral_collision_isomorphism_refinement"])
        self.assertFalse(self.acceptance["accepted_as_canonical_label_dataset"])
        self.assertFalse(self.acceptance["accepted_as_strict_spectral_source_law"])
        self.assertFalse(self.acceptance["accepted_as_ltotal_or_toe_promotion"])
        self.assertIn("canonical_label_dataset_exported", self.acceptance["missing_criteria"])
        self.assertFalse(self.payload["decision"]["negative_export_flags"]["toe_closure_exported"])

    def test_written_documents_reference_guardrail(self):
        self.assertIn("P2805/S1755", MD_PATH.read_text(encoding="utf-8"))
        self.assertIn("P2805/S1755", (ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md").read_text(encoding="utf-8"))
        self.assertIn("P2805/S1755", (ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md").read_text(encoding="utf-8"))
        self.assertIn("P2805", (REPO / "AGENTS.md").read_text(encoding="utf-8"))


if __name__ == "__main__":
    unittest.main()
