import json, subprocess, sys, unittest
from pathlib import Path
ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
SCRIPT = ROOT / "p3017_s1967_time_observable_action_eom_source_obstruction.py"
OUT = ROOT / "generated" / "p3017_s1967_time_observable_action_eom_source_obstruction.json"
MD = ROOT / "generated" / "p3017_s1967_time_observable_action_eom_source_obstruction.md"

class P3017TimeObservableActionEOMSourceTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        subprocess.run([sys.executable, str(SCRIPT)], cwd=ROOT, check=True, stdout=subprocess.PIPE, text=True)
        cls.payload = json.loads(OUT.read_text(encoding="utf-8"))

    def test_status_and_input(self):
        self.assertEqual(self.payload["status"], "P3017_TIME_OBSERVABLE_ACTION_EOM_SOURCE_OBSTRUCTION_NO_CLOSURE")
        self.assertIsNotNone(self.payload["input_hashes"]["P3016"])

    def test_formal_eom_and_unit_obstruction(self):
        cert = self.payload["finite_certificate"]
        self.assertEqual(cert["laplacian_rank"], 11)
        self.assertEqual(cert["laplacian_nullity"], 1)
        self.assertEqual(cert["residual_row_count"], 12)
        self.assertGreater(cert["nonzero_residual_count"], 0)
        self.assertEqual(cert["unit_obligation_count"], 4)
        self.assertEqual(cert["satisfied_unit_obligation_count"], 0)
        self.assertFalse(cert["accepted_as_unit_bearing_action_eom_source"])

    def test_obligations_and_negative_exports(self):
        obj = self.payload["constructed_theoretical_objects"]
        self.assertEqual(obj["object"], "TimeObservableQuadraticActionEOM_UnitSourceObstructionMatrix")
        obligations = {row["obligation"]: row["satisfied"] for row in obj["proof_obligations"]}
        self.assertTrue(obligations["typed_observable_input"])
        self.assertTrue(obligations["formal_variational_eom"])
        self.assertFalse(obligations["observable_solves_formal_eom"])
        self.assertFalse(obligations["unit_bearing_action_source"])
        self.assertFalse(obligations["hamiltonian_export"])
        for value in self.payload["decision"]["negative_export_flags"].values():
            self.assertFalse(value)

    def test_docs(self):
        self.assertIn("P3017/S1967", MD.read_text(encoding="utf-8"))
        self.assertIn("P3017/S1967", (ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md").read_text(encoding="utf-8"))
        self.assertIn("P3017/S1967", (ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md").read_text(encoding="utf-8"))
        self.assertIn("P3017", (REPO / "AGENTS.md").read_text(encoding="utf-8"))

if __name__ == "__main__":
    unittest.main()
