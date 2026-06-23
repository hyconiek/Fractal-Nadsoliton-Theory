import json, subprocess, sys, unittest
from pathlib import Path
ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
SCRIPT = ROOT / "p3055_s2005_mechanism_agnostic_selector_clue_sheaf.py"
OUT = ROOT / "generated" / "p3055_s2005_mechanism_agnostic_selector_clue_sheaf.json"
MD = ROOT / "generated" / "p3055_s2005_mechanism_agnostic_selector_clue_sheaf.md"

class P3055MechanismAgnosticSelectorClueSheafTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        subprocess.run([sys.executable, str(SCRIPT)], cwd=ROOT, check=True, stdout=subprocess.PIPE, text=True)
        cls.payload = json.loads(OUT.read_text(encoding="utf-8"))

    def test_status_and_input(self):
        self.assertEqual(self.payload["status"], "P3055_MECHANISM_AGNOSTIC_SELECTOR_CLUE_SHEAF_BOUNDARY_NO_EXPORT")
        self.assertIsNotNone(self.payload["input_hashes"]["P3054"])

    def test_finite_certificate(self):
        cert = self.payload["finite_certificate"]
        self.assertEqual(cert["content_grep_lanes"], 5)
        self.assertGreater(cert["content_grep_hits"], 0)
        self.assertEqual(cert["clue_stalks"], 5)
        self.assertEqual(cert["positive_clues"], 13)
        self.assertEqual(cert["obligation_atoms"], 6)
        self.assertEqual(cert["minimal_full_gluing_bundle_size"], 6)
        self.assertEqual(cert["minimal_full_gluing_bundle_count"], 1)
        self.assertEqual(cert["currently_satisfied_gluing_obligations"], 1)
        self.assertEqual(cert["proof_obligations"], 5)
        self.assertEqual(cert["satisfied_proof_obligations"], 3)
        self.assertFalse(cert["unknown_mechanism_selector_source_exported"])

    def test_sheaf_and_negative_exports(self):
        obj = self.payload["constructed_theoretical_objects"]
        self.assertEqual(obj["object"], "MechanismAgnosticSelectorClueSheaf")
        self.assertEqual(len(obj["stalks"]), 5)
        self.assertEqual(len(obj["obligations"]), 6)
        self.assertEqual(len(obj["minimal_full_gluing_bundles"]), 1)
        self.assertFalse(obj["unknown_mechanism_gluing_theorem_exported"])
        self.assertTrue(obj["current_artifact_satisfaction"]["nonzero_signed_value"])
        self.assertFalse(obj["current_artifact_satisfaction"]["strict_source_law"])
        for value in self.payload["decision"]["negative_export_flags"].values():
            self.assertFalse(value)

    def test_docs(self):
        self.assertIn("P3055/S2005", MD.read_text(encoding="utf-8"))
        self.assertIn("P3055/S2005", (ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md").read_text(encoding="utf-8"))
        self.assertIn("P3055/S2005", (ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md").read_text(encoding="utf-8"))
        self.assertIn("P3055", (REPO / "AGENTS.md").read_text(encoding="utf-8"))

if __name__ == "__main__":
    unittest.main()
