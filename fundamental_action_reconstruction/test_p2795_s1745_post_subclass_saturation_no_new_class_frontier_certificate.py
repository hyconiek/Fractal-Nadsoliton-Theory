import json
import unittest
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]
ROOT = REPO / "fundamental_action_reconstruction"
GEN = ROOT / "generated"
JSON_PATH = GEN / "p2795_s1745_post_subclass_saturation_no_new_class_frontier_certificate.json"
MD_PATH = GEN / "p2795_s1745_post_subclass_saturation_no_new_class_frontier_certificate.md"


class P2795PostSubclassSaturationNoNewClassFrontierCertificateTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        if not JSON_PATH.exists():
            import subprocess
            subprocess.run(["python", str(ROOT / "p2795_s1745_post_subclass_saturation_no_new_class_frontier_certificate.py")], check=True, cwd=ROOT)
        cls.payload = json.loads(JSON_PATH.read_text(encoding="utf-8"))
        cls.witness = cls.payload["post_subclass_saturation_witness"]
        cls.acceptance = cls.payload["acceptance_matrix"]

    def test_status_and_inputs(self):
        self.assertEqual(self.payload["status"], "P2795_POST_SUBCLASS_SATURATION_NO_NEW_CLASS_FRONTIER_CERTIFICATE_NO_CLOSURE")
        self.assertEqual(self.payload["input_statuses"]["P2791"], "P2791_EIGHT_CLASS_ORBIT_LOWER_BOUND_CERTIFICATE_NO_CLOSURE")
        self.assertEqual(self.payload["input_statuses"]["P2792"], "P2792_C16_TWO_JUMP_CIRCULANT_SUBCLASS_EXHAUSTION_CERTIFICATE_NO_CLOSURE")
        self.assertEqual(self.payload["input_statuses"]["P2793"], "P2793_Z2_4_FOUR_GENERATOR_CAYLEY_SUBCLASS_EXHAUSTION_CERTIFICATE_NO_CLOSURE")
        self.assertEqual(self.payload["input_statuses"]["P2794"], "P2794_Z2_4_GL4_TRANSITIVITY_AUTOMORPHISM_FACTORIZATION_CERTIFICATE_NO_CLOSURE")

    def test_named_subclasses_add_zero_new_classes_beyond_p2791(self):
        self.assertEqual(self.witness["p2791_eight_class_count"], 8)
        self.assertEqual(self.witness["union_named_subclass_label_count"], 6)
        self.assertEqual(self.witness["union_of_named_subclass_labels"], [
            "circulant_pm1_pm2",
            "circulant_pm1_pm3",
            "circulant_pm1_pm4",
            "circulant_pm1_pm6",
            "circulant_pm1_pm7",
            "torus_4x4",
        ])
        self.assertEqual(self.witness["uncovered_p2791_labels_by_named_subclasses"], ["p2790_eighth_witness", "two_c8_layers_cross_pm0_pm4"])
        self.assertEqual(self.witness["uncovered_p2791_label_count"], 2)
        self.assertTrue(self.witness["all_named_subclass_labels_inside_p2791"])
        self.assertEqual(self.witness["total_new_labels_beyond_p2791"], 0)
        self.assertEqual(self.witness["total_new_orbit_lower_bound_added_beyond_p2791"], 0)

    def test_subclass_rows_are_saturated(self):
        rows = self.witness["subclass_rows"]
        self.assertEqual([row["certified_class_count"] for row in rows], [5, 1, 1])
        self.assertTrue(all(row["all_matches_inside_p2791_eight_class_set"] for row in rows))
        self.assertTrue(all(row["new_labels_not_already_in_p2791"] == [] for row in rows))
        self.assertEqual([row["marginal_new_p2791_label_count"] for row in rows], [5, 1, 0])

    def test_acceptance_blocks_closure(self):
        self.assertTrue(self.acceptance["accepted_as_post_subclass_saturation_no_new_class_certificate"])
        self.assertFalse(self.acceptance["accepted_as_full_16node_canonical_generator_certificate"])
        self.assertFalse(self.acceptance["accepted_as_strict_spectral_source_law"])
        self.assertFalse(self.acceptance["accepted_as_canonical_nadsoliton_geometry_source"])
        self.assertFalse(self.acceptance["accepted_as_ltotal_or_toe_promotion"])
        self.assertIn("P2697-P2795", self.payload["decision"]["next_honest_step"])
        self.assertTrue(all(value is False for value in self.payload["decision"]["negative_export_flags"].values()))

    def test_documentation_updated(self):
        self.assertIn("P2795/S1745", MD_PATH.read_text(encoding="utf-8"))
        self.assertIn("P2795/S1745", (ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md").read_text(encoding="utf-8"))
        self.assertIn("P2795/S1745", (ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md").read_text(encoding="utf-8"))
        self.assertIn("P2795", (REPO / "AGENTS.md").read_text(encoding="utf-8"))


if __name__ == "__main__":
    unittest.main()
