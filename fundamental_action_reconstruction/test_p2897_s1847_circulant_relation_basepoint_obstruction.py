import json
import subprocess
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
SCRIPT = ROOT / "p2897_s1847_circulant_relation_basepoint_obstruction.py"
JSON_PATH = ROOT / "generated" / "p2897_s1847_circulant_relation_basepoint_obstruction.json"
MD_PATH = ROOT / "generated" / "p2897_s1847_circulant_relation_basepoint_obstruction.md"


class P2897CirculantRelationBasepointObstructionTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        subprocess.run([sys.executable, str(SCRIPT)], cwd=ROOT, check=True, stdout=subprocess.DEVNULL)
        cls.payload = json.loads(JSON_PATH.read_text(encoding="utf-8"))
        cls.audit = cls.payload["circulant_relation_basepoint_obstruction"]["relation_exhaustion"]
        cls.facts = cls.payload["acceptance_matrix"]["facts"]
        cls.flags = cls.payload["decision"]["negative_export_flags"]

    def test_status_and_input(self):
        self.assertEqual(self.payload["status"], "P2897_CIRCULANT_RELATION_BASEPOINT_OBSTRUCTION_NO_CLOSURE")
        self.assertTrue(self.facts["p2896_rechecked"])
        self.assertEqual(
            self.payload["circulant_relation_basepoint_obstruction"]["input_status_rechecked"],
            "P2896_INVARIANT_SCALAR_BASEPOINT_LAW_EXHAUSTION_NO_CLOSURE",
        )

    def test_relation_exhaustion_counts(self):
        self.assertEqual(self.audit["torsor_size"], 12)
        self.assertEqual(self.audit["directed_circulant_relation_count"], 4096)
        self.assertEqual(self.audit["undirected_loopless_circulant_graph_count"], 64)
        self.assertTrue(self.audit["all_relations_have_translation_automorphisms"])
        self.assertEqual(self.audit["translation_orbit_size_histogram"], {"12": 4096})
        self.assertEqual(self.audit["relations_selecting_unique_vertex_count"], 0)
        self.assertEqual(self.audit["relations_with_unique_degree_vertex_count"], 0)
        self.assertEqual(self.audit["relations_with_unique_local_row_profile_count"], 0)

    def test_false_export_flags(self):
        self.assertFalse(any(self.flags.values()))
        self.assertFalse(self.payload["acceptance_matrix"]["accepted_as_nonimported_basepoint_law"])
        self.assertFalse(self.payload["acceptance_matrix"]["exports_translation_breaking_relation_source"])
        self.assertFalse(self.payload["acceptance_matrix"]["exports_unique_coupling_to_9_over_5_carrier"])

    def test_documents_updated(self):
        self.assertIn("P2897/S1847", MD_PATH.read_text(encoding="utf-8"))
        self.assertIn("P2897/S1847", (ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md").read_text(encoding="utf-8"))
        self.assertIn("P2897/S1847", (ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md").read_text(encoding="utf-8"))
        self.assertIn("P2897", (REPO / "AGENTS.md").read_text(encoding="utf-8"))


if __name__ == "__main__":
    unittest.main()
