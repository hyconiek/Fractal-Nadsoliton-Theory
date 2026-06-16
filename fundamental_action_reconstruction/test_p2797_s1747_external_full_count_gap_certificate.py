import json
import unittest
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]
ROOT = REPO / "fundamental_action_reconstruction"
GEN = ROOT / "generated"
JSON_PATH = GEN / "p2797_s1747_external_full_count_gap_certificate.json"
MD_PATH = GEN / "p2797_s1747_external_full_count_gap_certificate.md"


class P2797ExternalFullCountGapCertificateTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        if not JSON_PATH.exists():
            import subprocess
            subprocess.run(["python", str(ROOT / "p2797_s1747_external_full_count_gap_certificate.py")], check=True, cwd=ROOT)
        cls.payload = json.loads(JSON_PATH.read_text(encoding="utf-8"))
        cls.witness = cls.payload["external_full_count_gap_witness"]
        cls.acceptance = cls.payload["acceptance_matrix"]

    def test_status_and_inputs(self):
        self.assertEqual(self.payload["status"], "P2797_EXTERNAL_FULL_COUNT_GAP_CERTIFICATE_NO_CLOSURE")
        self.assertEqual(self.payload["input_statuses"]["P2791"], "P2791_EIGHT_CLASS_ORBIT_LOWER_BOUND_CERTIFICATE_NO_CLOSURE")
        self.assertEqual(self.payload["input_statuses"]["P2795"], "P2795_POST_SUBCLASS_SATURATION_NO_NEW_CLASS_FRONTIER_CERTIFICATE_NO_CLOSURE")
        self.assertEqual(self.payload["input_statuses"]["P2796"], "P2796_FULL_GENERATOR_ARTIFACT_ABSENCE_PROVENANCE_GATE_NO_CLOSURE")

    def test_external_count_gap_numbers(self):
        self.assertEqual(self.witness["external_full_connected_16node_4regular_class_count"], 8037418)
        self.assertEqual(self.witness["p2791_distinct_class_lower_bound"], 8)
        self.assertEqual(self.witness["p2795_named_subclass_union_class_count"], 6)
        self.assertEqual(self.witness["p2795_new_classes_beyond_p2791"], 0)
        self.assertEqual(self.witness["current_distinct_class_lower_bound"], 8)
        self.assertEqual(self.witness["uncovered_class_gap_after_p2791"], 8037410)
        self.assertEqual(self.witness["uncovered_class_gap_after_named_subclasses"], 8037410)
        self.assertEqual(self.witness["p2791_coverage_fraction_exact"], "4/4018709")
        self.assertEqual(self.witness["named_subclass_union_fraction_exact"], "3/4018709")

    def test_external_source_is_scoped(self):
        source = self.witness["external_source"]
        self.assertEqual(source["url"], "https://www.mathe2.uni-bayreuth.de/markus/reggraphs.html")
        self.assertEqual(source["table"], "Connected regular graphs")
        self.assertEqual(source["row"], "vertices=16")
        self.assertEqual(source["column"], "degree=4")
        self.assertEqual(source["reported_count"], 8037418)

    def test_acceptance_blocks_closure(self):
        self.assertTrue(self.acceptance["accepted_as_external_full_count_gap_certificate"])
        self.assertFalse(self.acceptance["accepted_as_full_16node_canonical_generator_certificate"])
        self.assertFalse(self.acceptance["accepted_as_strict_spectral_source_law"])
        self.assertFalse(self.acceptance["accepted_as_canonical_nadsoliton_geometry_source"])
        self.assertFalse(self.acceptance["accepted_as_ltotal_or_toe_promotion"])
        self.assertIn("P2697-P2797", self.payload["decision"]["next_honest_step"])
        self.assertTrue(all(value is False for value in self.payload["decision"]["negative_export_flags"].values()))

    def test_documentation_updated(self):
        self.assertIn("P2797/S1747", MD_PATH.read_text(encoding="utf-8"))
        self.assertIn("P2797/S1747", (ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md").read_text(encoding="utf-8"))
        self.assertIn("P2797/S1747", (ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md").read_text(encoding="utf-8"))
        self.assertIn("P2797", (REPO / "AGENTS.md").read_text(encoding="utf-8"))


if __name__ == "__main__":
    unittest.main()
