import json
import subprocess
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
SCRIPT = ROOT / "p2899_s1849_post_defect_potential_readiness_matrix.py"
JSON_PATH = ROOT / "generated" / "p2899_s1849_post_defect_potential_readiness_matrix.json"
MD_PATH = ROOT / "generated" / "p2899_s1849_post_defect_potential_readiness_matrix.md"


class P2899PostDefectPotentialReadinessMatrixTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        subprocess.run([sys.executable, str(SCRIPT)], cwd=ROOT, check=True, stdout=subprocess.DEVNULL)
        cls.payload = json.loads(JSON_PATH.read_text(encoding="utf-8"))
        cls.matrix = cls.payload["potential_readiness_matrix"]
        cls.acceptance = cls.payload["acceptance_matrix"]
        cls.flags = cls.payload["decision"]["negative_export_flags"]

    def test_status_and_inputs(self):
        self.assertEqual(self.payload["status"], "P2899_POST_DEFECT_POTENTIAL_READINESS_MATRIX_NO_CLOSURE")
        self.assertTrue(self.acceptance["all_inputs_rechecked"])
        self.assertEqual(set(self.payload["input_statuses"]), {"P2895", "P2896", "P2897", "P2898"})

    def test_potential_and_blocker_counts(self):
        self.assertEqual(self.matrix["positive_symptom_count"], 4)
        self.assertEqual(self.matrix["blocker_count"], 6)
        self.assertFalse(self.matrix["closure_ready"])
        self.assertEqual(self.matrix["toe_potential_class"], "conditional_structural_potential_but_no_current_ToE_closure")
        self.assertTrue(all(row["positive"] for row in self.matrix["positive_symptoms"]))
        self.assertTrue(all(row["blocked"] for row in self.matrix["blockers"]))

    def test_acceptance_boundary(self):
        self.assertTrue(self.acceptance["positive_symptoms_exist"])
        self.assertTrue(self.acceptance["accepted_as_toe_potential_evidence"])
        self.assertFalse(self.acceptance["accepted_as_toe_closure"])
        self.assertFalse(any(self.flags.values()))

    def test_documents_updated(self):
        self.assertIn("P2899/S1849", MD_PATH.read_text(encoding="utf-8"))
        self.assertIn("P2899/S1849", (ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md").read_text(encoding="utf-8"))
        self.assertIn("P2899/S1849", (ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md").read_text(encoding="utf-8"))
        self.assertIn("P2899", (REPO / "AGENTS.md").read_text(encoding="utf-8"))


if __name__ == "__main__":
    unittest.main()
