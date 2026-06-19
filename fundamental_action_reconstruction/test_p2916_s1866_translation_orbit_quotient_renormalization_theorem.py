import json
import subprocess
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
SCRIPT = ROOT / "p2916_s1866_translation_orbit_quotient_renormalization_theorem.py"
JSON_PATH = ROOT / "generated" / "p2916_s1866_translation_orbit_quotient_renormalization_theorem.json"
MD_PATH = ROOT / "generated" / "p2916_s1866_translation_orbit_quotient_renormalization_theorem.md"


class P2916TranslationOrbitQuotientRenormalizationTheoremTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        subprocess.run([sys.executable, str(SCRIPT)], cwd=ROOT, check=True, stdout=subprocess.DEVNULL)
        cls.payload = json.loads(JSON_PATH.read_text(encoding="utf-8"))
        cls.acceptance = cls.payload["acceptance_matrix"]
        cls.objects = cls.payload["constructed_theoretical_objects"]
        cls.flags = cls.payload["decision"]["negative_export_flags"]

    def test_status_and_p2915_input(self):
        self.assertEqual(self.payload["status"], "P2916_TRANSLATION_ORBIT_QUOTIENT_RENORMALIZATION_THEOREM_FINITE_EXPORT_NO_LTOTAL")
        self.assertTrue(self.acceptance["p2915_rechecked_three_candidate_boundary"])

    def test_translation_orbit_quotient_selected(self):
        self.assertEqual(self.acceptance["selected_quotient"], "displacement_quotient_q(i,j)=j-i_mod12")
        self.assertEqual(self.acceptance["translation_orbit_count"], 12)
        self.assertTrue(self.acceptance["all_translation_orbits_size_12"])
        self.assertEqual(self.acceptance["displacement_label_translation_failure_count"], 0)
        self.assertGreater(self.acceptance["source_label_translation_failure_count"], 0)
        self.assertGreater(self.acceptance["target_label_translation_failure_count"], 0)

    def test_renormalization_and_lay_explanation(self):
        self.assertEqual(self.acceptance["per_edge_renormalized_weight"], "1/144")
        self.assertEqual(self.acceptance["quotient_total_weight"], "1/1")
        self.assertTrue(self.acceptance["finite_translation_quotient_theorem_exported"])
        self.assertIn("12 Z12", self.payload["lay_explanation"]["what_12_counts"])
        self.assertIn("12 x 12", self.payload["lay_explanation"]["what_144_counts"])

    def test_no_ltotal_or_closure_export(self):
        self.assertFalse(self.acceptance["strict_gamma_9_5_source_exported"])
        self.assertFalse(self.acceptance["continuum_field_variable_provenance_exported"])
        self.assertFalse(self.acceptance["accepted_as_nonproxy_ltotal_measure_bridge"])
        self.assertFalse(any(self.flags.values()))

    def test_documents_updated(self):
        self.assertIn("P2916/S1866", MD_PATH.read_text(encoding="utf-8"))
        self.assertIn("P2916/S1866", (ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md").read_text(encoding="utf-8"))
        self.assertIn("P2916/S1866", (ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md").read_text(encoding="utf-8"))
        self.assertIn("P2916", (REPO / "AGENTS.md").read_text(encoding="utf-8"))


if __name__ == "__main__":
    unittest.main()
