import json, subprocess, sys, unittest
from pathlib import Path
ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
SCRIPT = ROOT / "p3018_s1968_lambda_action_unit_normalization_candidate_obstruction.py"
OUT = ROOT / "generated" / "p3018_s1968_lambda_action_unit_normalization_candidate_obstruction.json"
MD = ROOT / "generated" / "p3018_s1968_lambda_action_unit_normalization_candidate_obstruction.md"

class P3018LambdaActionUnitNormalizationTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        subprocess.run([sys.executable, str(SCRIPT)], cwd=ROOT, check=True, stdout=subprocess.PIPE, text=True)
        cls.payload = json.loads(OUT.read_text(encoding="utf-8"))

    def test_status_and_input(self):
        self.assertEqual(self.payload["status"], "P3018_LAMBDA_ACTION_UNIT_NORMALIZATION_CANDIDATE_OBSTRUCTION_NO_CLOSURE")
        self.assertIsNotNone(self.payload["input_hashes"]["P3017"])

    def test_positive_candidate_but_scale_obstruction(self):
        cert = self.payload["finite_certificate"]
        self.assertTrue(cert["gradient_energy_positive"])
        self.assertTrue(cert["lambda_star_positive"])
        self.assertAlmostEqual(cert["normalized_action"], 1.0, places=12)
        self.assertEqual(cert["scale_row_count"], 4)
        self.assertEqual(cert["scale_invariant_lambda_rows"], 1)
        self.assertFalse(cert["accepted_as_strict_lambda_action_unit_source"])

    def test_obligations_and_negative_exports(self):
        obj = self.payload["constructed_theoretical_objects"]
        self.assertEqual(obj["object"], "LambdaActionUnitNormalization_CandidateObstructionMatrix")
        obligations = {row["obligation"]: row["satisfied"] for row in obj["proof_obligations"]}
        self.assertTrue(obligations["positive_lambda_candidate"])
        self.assertTrue(obligations["normalizes_formal_action_to_one"])
        self.assertFalse(obligations["observable_scale_invariant_lambda"])
        self.assertFalse(obligations["strict_action_quantum_source"])
        self.assertFalse(obligations["strict_observable_unit_source"])
        self.assertFalse(obligations["clock_and_hamiltonian_unit_source"])
        for value in self.payload["decision"]["negative_export_flags"].values():
            self.assertFalse(value)

    def test_docs(self):
        self.assertIn("P3018/S1968", MD.read_text(encoding="utf-8"))
        self.assertIn("P3018/S1968", (ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md").read_text(encoding="utf-8"))
        self.assertIn("P3018/S1968", (ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md").read_text(encoding="utf-8"))
        self.assertIn("P3018", (REPO / "AGENTS.md").read_text(encoding="utf-8"))

if __name__ == "__main__":
    unittest.main()
