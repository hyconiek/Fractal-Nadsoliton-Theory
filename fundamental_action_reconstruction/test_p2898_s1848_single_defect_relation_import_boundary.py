import json
import subprocess
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
SCRIPT = ROOT / "p2898_s1848_single_defect_relation_import_boundary.py"
JSON_PATH = ROOT / "generated" / "p2898_s1848_single_defect_relation_import_boundary.json"
MD_PATH = ROOT / "generated" / "p2898_s1848_single_defect_relation_import_boundary.md"


class P2898SingleDefectRelationImportBoundaryTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        subprocess.run([sys.executable, str(SCRIPT)], cwd=ROOT, check=True, stdout=subprocess.DEVNULL)
        cls.payload = json.loads(JSON_PATH.read_text(encoding="utf-8"))
        cls.audit = cls.payload["single_defect_relation_import_boundary"]["single_defect_audit"]
        cls.facts = cls.payload["acceptance_matrix"]["facts"]
        cls.flags = cls.payload["decision"]["negative_export_flags"]

    def test_status_and_input(self):
        self.assertEqual(self.payload["status"], "P2898_SINGLE_DEFECT_RELATION_IMPORT_BOUNDARY_NO_CLOSURE")
        self.assertTrue(self.facts["p2897_rechecked"])
        self.assertEqual(
            self.payload["single_defect_relation_import_boundary"]["input_status_rechecked"],
            "P2897_CIRCULANT_RELATION_BASEPOINT_OBSTRUCTION_NO_CLOSURE",
        )

    def test_single_defect_counts(self):
        self.assertEqual(self.audit["torsor_size"], 12)
        self.assertEqual(self.audit["circulant_background_count"], 4096)
        self.assertEqual(self.audit["directed_edge_count"], 144)
        self.assertEqual(self.audit["single_defect_candidate_count"], 589824)
        self.assertEqual(self.audit["translation_quotient_orbit_count"], 49152)
        self.assertEqual(self.audit["defect_orbit_size_histogram"], {"12": 49152})
        self.assertEqual(self.audit["source_neutral_defect_placement_count"], 0)
        self.assertTrue(self.audit["all_defect_orbits_free"])

    def test_false_export_flags(self):
        self.assertFalse(any(self.flags.values()))
        self.assertFalse(self.payload["acceptance_matrix"]["accepted_as_nonimported_basepoint_law"])
        self.assertFalse(self.payload["acceptance_matrix"]["exports_translation_breaking_relation_source"])
        self.assertFalse(self.payload["acceptance_matrix"]["exports_unique_coupling_to_9_over_5_carrier"])

    def test_documents_updated(self):
        self.assertIn("P2898/S1848", MD_PATH.read_text(encoding="utf-8"))
        self.assertIn("P2898/S1848", (ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md").read_text(encoding="utf-8"))
        self.assertIn("P2898/S1848", (ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md").read_text(encoding="utf-8"))
        self.assertIn("P2898", (REPO / "AGENTS.md").read_text(encoding="utf-8"))


if __name__ == "__main__":
    unittest.main()
