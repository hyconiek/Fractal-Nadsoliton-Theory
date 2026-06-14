import json
import subprocess
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
SCRIPT = ROOT / "p2737_s1687_lay_toe_potential_readiness_matrix.py"
OUT = ROOT / "generated" / "p2737_s1687_lay_toe_potential_readiness_matrix.json"
MD = ROOT / "generated" / "p2737_s1687_lay_toe_potential_readiness_matrix.md"


class P2737LayToePotentialReadinessMatrixTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        subprocess.run([sys.executable, str(SCRIPT)], cwd=ROOT, check=True)
        cls.payload = json.loads(OUT.read_text(encoding="utf-8"))
        cls.scan = cls.payload["content_evidence_scan"]
        cls.matrix = cls.payload["toe_readiness_matrix"]
        cls.acceptance = cls.payload["acceptance_matrix"]

    def test_content_scan_has_hits_for_all_toe_lanes(self):
        self.assertEqual(self.payload["status"], "P2737_LAY_TOE_POTENTIAL_READINESS_MATRIX_NO_CLOSURE")
        self.assertEqual(self.scan["content_pattern_count"], 7)
        self.assertTrue(self.scan["all_patterns_have_hits"])
        self.assertGreater(self.scan["hit_counts"]["toe_blocker"], 0)
        self.assertGreater(self.scan["hit_counts"]["selector_blocker"], 0)

    def test_readiness_fraction_supports_potential_not_closure(self):
        self.assertEqual(self.matrix["readiness_fraction"], "2/8")
        self.assertFalse(self.matrix["closure_allowed"])
        met = {row["condition"]: row["met"] for row in self.matrix["rows"]}
        self.assertTrue(met["corrected_ontology"])
        self.assertTrue(met["finite_computational_discipline"])
        self.assertFalse(met["strict_selector_or_orientation_source"])
        self.assertFalse(met["toe_closure"])

    def test_acceptance_blocks_toe_closure(self):
        self.assertFalse(self.acceptance["accepted_as_toe_closure"])
        self.assertIn("strict_selector_source_closed", self.acceptance["missing_criteria"])
        self.assertIn("bridge_completion_closed", self.acceptance["missing_criteria"])
        self.assertIn("toe_closure_allowed", self.acceptance["missing_criteria"])
        self.assertTrue(all(value is False for value in self.payload["decision"]["negative_export_flags"].values()))

    def test_lay_explanation_and_recommendation_are_plain_and_guarded(self):
        explanation = self.payload["lay_explanation"]
        recommendation = self.payload["decision"]["next_honest_step"]
        self.assertIn("compass needle", explanation)
        self.assertIn("not a ToE yet", explanation)
        self.assertIn("do not make a closure claim", recommendation)
        self.assertIn("P2697-P2737", recommendation)

    def test_documentation_updated(self):
        self.assertIn("P2737/S1687", MD.read_text(encoding="utf-8"))
        self.assertIn("P2737/S1687", (ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md").read_text(encoding="utf-8"))
        self.assertIn("P2737/S1687", (ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md").read_text(encoding="utf-8"))
        self.assertIn("P2737", (REPO / "AGENTS.md").read_text(encoding="utf-8"))


if __name__ == "__main__":
    unittest.main()
