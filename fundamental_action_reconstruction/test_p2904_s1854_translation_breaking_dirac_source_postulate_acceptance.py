import json
import subprocess
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
SCRIPT = ROOT / "p2904_s1854_translation_breaking_dirac_source_postulate_acceptance.py"
JSON_PATH = ROOT / "generated" / "p2904_s1854_translation_breaking_dirac_source_postulate_acceptance.json"
MD_PATH = ROOT / "generated" / "p2904_s1854_translation_breaking_dirac_source_postulate_acceptance.md"


class P2904TranslationBreakingDiracSourcePostulateAcceptanceTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        subprocess.run([sys.executable, str(SCRIPT)], cwd=ROOT, check=True, stdout=subprocess.DEVNULL)
        cls.payload = json.loads(JSON_PATH.read_text(encoding="utf-8"))
        cls.acceptance = cls.payload["acceptance_matrix"]
        cls.objects = cls.payload["constructed_theoretical_objects"]
        cls.flags = cls.payload["decision"]["negative_export_flags"]

    def test_status_and_input(self):
        self.assertEqual(self.payload["status"], "P2904_TRANSLATION_BREAKING_DIRAC_SOURCE_POSTULATE_ACCEPTED_AS_CANDIDATE_NO_STRICT_PROVENANCE")
        self.assertTrue(self.acceptance["p2903_rechecked"])
        self.assertTrue(self.acceptance["translation_breaking_source_constructed"])
        self.assertTrue(self.acceptance["passes_p2903_fixed_point_obstruction"])

    def test_source_is_translation_breaking_and_coupled(self):
        self.assertEqual(self.objects["source_values"], [1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0])
        self.assertEqual(self.acceptance["nonzero_computed_value"], 1)
        self.assertEqual(self.acceptance["unique_signed_support_count"], 1)
        self.assertEqual(self.acceptance["translation_stabilizer_size"], 1)
        self.assertEqual(self.acceptance["translation_orbit_size"], 12)
        self.assertTrue(self.acceptance["coupling_to_one_pointed_axiom_constructed"])
        self.assertEqual(self.acceptance["selected_basepoint"], 0)
        self.assertEqual(self.acceptance["selected_polarity"], 1)
        self.assertEqual(self.acceptance["defect_edge"], [0, 5])

    def test_remaining_strict_blocks(self):
        self.assertFalse(self.acceptance["strict_nadsoliton_provenance_exported"])
        self.assertFalse(self.acceptance["unit_bearing_nonproxy_ltotal_coupling_exported"])
        self.assertFalse(any(self.flags.values()))
        self.assertFalse(self.flags["nonproxy_ltotal_exported"])
        self.assertFalse(self.flags["toe_closure_exported"])

    def test_documents_updated(self):
        self.assertIn("P2904/S1854", MD_PATH.read_text(encoding="utf-8"))
        self.assertIn("P2904/S1854", (ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md").read_text(encoding="utf-8"))
        self.assertIn("P2904/S1854", (ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md").read_text(encoding="utf-8"))
        self.assertIn("P2904", (REPO / "AGENTS.md").read_text(encoding="utf-8"))


if __name__ == "__main__":
    unittest.main()
