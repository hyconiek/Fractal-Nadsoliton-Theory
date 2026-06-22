import json, subprocess, sys, unittest
from pathlib import Path
ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
SCRIPT = ROOT / "p3011_s1961_toe_strict_kernel_selector_role_separation_matrix.py"
OUT = ROOT / "generated" / "p3011_s1961_toe_strict_kernel_selector_role_separation_matrix.json"
MD = ROOT / "generated" / "p3011_s1961_toe_strict_kernel_selector_role_separation_matrix.md"

class P3011ToEStrictKernelSelectorRoleSeparationTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        subprocess.run([sys.executable, str(SCRIPT)], cwd=ROOT, check=True, stdout=subprocess.PIPE, text=True)
        cls.payload = json.loads(OUT.read_text(encoding="utf-8"))

    def test_status_and_input(self):
        self.assertEqual(self.payload["status"], "P3011_TOE_STRICT_KERNEL_SELECTOR_ROLE_SEPARATION_MATRIX_NO_CLOSURE")
        self.assertIsNotNone(self.payload["input_hashes"]["P3010"])

    def test_kernel_and_selector_certificate(self):
        cert = self.payload["answer_certificate"]
        self.assertIn("operational gate-selected correlation/compression profile over separation d", cert["strict_kernel_describes"])
        self.assertIn("directed representative", cert["selector_is"])
        self.assertEqual(cert["orientation_reversing_units_on_pair"], [11])
        self.assertEqual(cert["aut_invariant_directed_choice_count"], 0)
        self.assertEqual(cert["acceptance_matrix_rows"], 128)
        self.assertEqual(cert["accepted_rows"], 1)

    def test_role_rows_and_nonpromotion(self):
        obj = self.payload["constructed_theoretical_objects"]
        self.assertEqual(obj["object"], "ToEStrictKernelSelectorRoleSeparation_ObligationMatrix")
        rows = {r["row"]: r for r in obj["role_separation_rows"]}
        self.assertTrue(rows["ontology"]["satisfied"])
        self.assertTrue(rows["strict_kernel"]["satisfied"])
        self.assertFalse(rows["selector"]["satisfied"])
        self.assertFalse(rows["toe_closure"]["satisfied"])
        for value in self.payload["decision"]["negative_export_flags"].values():
            self.assertFalse(value)

    def test_docs(self):
        self.assertIn("P3011/S1961", MD.read_text(encoding="utf-8"))
        self.assertIn("P3011/S1961", (ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md").read_text(encoding="utf-8"))
        self.assertIn("P3011/S1961", (ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md").read_text(encoding="utf-8"))
        self.assertIn("P3011", (REPO / "AGENTS.md").read_text(encoding="utf-8"))

if __name__ == "__main__":
    unittest.main()
