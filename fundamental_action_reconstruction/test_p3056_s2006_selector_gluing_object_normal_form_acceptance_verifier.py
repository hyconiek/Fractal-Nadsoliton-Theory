import json, subprocess, sys, unittest
from pathlib import Path
ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
SCRIPT = ROOT / "p3056_s2006_selector_gluing_object_normal_form_acceptance_verifier.py"
OUT = ROOT / "generated" / "p3056_s2006_selector_gluing_object_normal_form_acceptance_verifier.json"
MD = ROOT / "generated" / "p3056_s2006_selector_gluing_object_normal_form_acceptance_verifier.md"

class P3056SelectorGluingObjectNormalFormTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        subprocess.run([sys.executable, str(SCRIPT)], cwd=ROOT, check=True, stdout=subprocess.PIPE, text=True)
        cls.payload = json.loads(OUT.read_text(encoding="utf-8"))

    def test_status_and_input(self):
        self.assertEqual(self.payload["status"], "P3056_SELECTOR_GLUING_OBJECT_NORMAL_FORM_ACCEPTANCE_VERIFIER_NO_EXPORT")
        self.assertIsNotNone(self.payload["input_hashes"]["P3055"])

    def test_finite_certificate(self):
        cert = self.payload["finite_certificate"]
        self.assertEqual(cert["content_grep_lanes"], 3)
        self.assertGreater(cert["content_grep_hits"], 0)
        self.assertEqual(cert["typed_rows"], 6)
        self.assertEqual(cert["compatibility_squares"], 6)
        self.assertEqual(cert["candidate_carriers"], 6)
        self.assertEqual(cert["carrier_subsets_enumerated"], 63)
        self.assertEqual(cert["row_complete_bundles"], 0)
        self.assertEqual(cert["accepted_gluing_objects"], 0)
        self.assertEqual(cert["proof_obligations"], 5)
        self.assertEqual(cert["satisfied_proof_obligations"], 3)

    def test_normal_form_and_negative_exports(self):
        obj = self.payload["constructed_theoretical_objects"]
        self.assertEqual(obj["object"], "SelectorGluingObjectNormalForm")
        self.assertEqual(len(obj["typed_rows"]), 6)
        self.assertEqual(len(obj["compatibility_squares"]), 6)
        self.assertEqual(len(obj["candidate_carriers"]), 6)
        self.assertEqual(len(obj["carrier_subset_pushout"]), 63)
        for row in obj["carrier_subset_pushout"]:
            self.assertFalse(row["accepted_gluing_object"])
        for value in self.payload["decision"]["negative_export_flags"].values():
            self.assertFalse(value)

    def test_docs(self):
        self.assertIn("P3056/S2006", MD.read_text(encoding="utf-8"))
        self.assertIn("P3056/S2006", (ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md").read_text(encoding="utf-8"))
        self.assertIn("P3056/S2006", (ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md").read_text(encoding="utf-8"))
        self.assertIn("P3056", (REPO / "AGENTS.md").read_text(encoding="utf-8"))

if __name__ == "__main__":
    unittest.main()
