import json
import subprocess
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
SCRIPT = ROOT / "p2893_s1843_free_orbit_section_obstruction_theorem.py"
JSON_PATH = ROOT / "generated" / "p2893_s1843_free_orbit_section_obstruction_theorem.json"
MD_PATH = ROOT / "generated" / "p2893_s1843_free_orbit_section_obstruction_theorem.md"


class P2893FreeOrbitSectionObstructionTheoremTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        subprocess.run([sys.executable, str(SCRIPT)], cwd=ROOT, check=True, stdout=subprocess.DEVNULL)
        cls.payload = json.loads(JSON_PATH.read_text(encoding="utf-8"))
        cls.audit = cls.payload["free_orbit_section_obstruction_theorem"]["finite_section_audit"]
        cls.facts = cls.payload["acceptance_matrix"]["facts"]
        cls.flags = cls.payload["decision"]["negative_export_flags"]

    def test_status_and_input(self):
        self.assertEqual(self.payload["status"], "P2893_FREE_ORBIT_SECTION_OBSTRUCTION_THEOREM_NO_CLOSURE")
        self.assertTrue(self.facts["p2892_rechecked"])
        self.assertEqual(self.payload["free_orbit_section_obstruction_theorem"]["input_status_rechecked"], "P2892_POST_PHASE_ORIGIN_INVENTORY_STATE_MAP_NO_NEW_LIVE_FRONTIER_CERTIFICATE")

    def test_finite_section_counts(self):
        self.assertEqual(self.audit["target_triple_count"], 600)
        self.assertEqual(self.audit["translation_orbit_count"], 50)
        self.assertEqual(self.audit["orbit_size_histogram"], {"12": 50})
        self.assertEqual(self.audit["point_stabilizer_size_histogram_global"], {"1": 600})
        self.assertEqual(self.audit["quotient_to_embedded_invariant_section_count_total"], 0)
        self.assertEqual(self.audit["equivariant_endomorphism_count_total_not_source_neutral"], 600)

    def test_false_export_flags(self):
        self.assertFalse(any(self.flags.values()))
        self.assertFalse(self.payload["acceptance_matrix"]["accepted_as_source_neutral_embedded_representative_law"])
        self.assertFalse(self.payload["acceptance_matrix"]["exports_translation_breaking_source"])
        self.assertFalse(self.payload["acceptance_matrix"]["exports_nonimported_9_over_5_variational_density"])

    def test_documents_updated(self):
        self.assertIn("P2893/S1843", MD_PATH.read_text(encoding="utf-8"))
        self.assertIn("P2893/S1843", (ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md").read_text(encoding="utf-8"))
        self.assertIn("P2893/S1843", (ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md").read_text(encoding="utf-8"))
        self.assertIn("P2893", (REPO / "AGENTS.md").read_text(encoding="utf-8"))


if __name__ == "__main__":
    unittest.main()
