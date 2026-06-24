import json, subprocess, sys, unittest
from pathlib import Path
ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
SCRIPT = ROOT / "p3064_s2014_strict_polarity_source_theorem_sat.py"
OUT = ROOT / "generated" / "p3064_s2014_strict_polarity_source_theorem_sat.json"
MD = ROOT / "generated" / "p3064_s2014_strict_polarity_source_theorem_sat.md"

class P3064StrictPolaritySourceTheoremSatTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        subprocess.run([sys.executable, str(SCRIPT)], cwd=ROOT, check=True, stdout=subprocess.PIPE, text=True)
        cls.payload = json.loads(OUT.read_text(encoding="utf-8"))

    def test_status_and_input(self):
        self.assertEqual(self.payload["status"], "P3064_STRICT_POLARITY_SOURCE_THEOREM_SAT_NO_CURRENT_EXPORT")
        self.assertIsNotNone(self.payload["input_hashes"]["P3063"])

    def test_content_first_and_sat_counts(self):
        cert = self.payload["finite_certificate"]
        self.assertEqual(cert["content_grep_lanes"], 3)
        self.assertGreater(cert["content_grep_hits"], 0)
        self.assertEqual(cert["primitive_atoms"], 4)
        self.assertEqual(cert["sat_rows"], 16)
        self.assertEqual(cert["consistent_rows"], 9)
        self.assertEqual(cert["accepted_theorem_rows"], 4)
        self.assertEqual(cert["minimal_accepting_atom_count"], 2)
        self.assertEqual(cert["minimal_accepting_rows"], 4)
        self.assertEqual(cert["current_exported_atoms"], 0)
        self.assertFalse(cert["current_row_accepted"])
        self.assertEqual(cert["satisfied_proof_obligations"], 3)

    def test_theorem_object_and_current_row(self):
        obj = self.payload["constructed_theoretical_objects"]
        theorem = obj["minimal_strict_polarity_source_theorem"]
        self.assertEqual(theorem["object"], "MinimalStrictPolaritySourceTheorem")
        self.assertEqual(len(theorem["primitive_atoms"]), 4)
        current = obj["current_artifact_row"]
        self.assertFalse(current["accepted_minimal_strict_polarity_source_theorem"])
        self.assertIn("absolute_nonpremise_source_sign", current["missing_or_blocking"])
        self.assertIn("absolute_nonpremise_coupling_polarity", current["missing_or_blocking"])
        self.assertEqual(len(obj["minimal_accepting_rows"]), 4)

    def test_negative_exports_and_docs(self):
        for value in self.payload["decision"]["negative_export_flags"].values():
            self.assertFalse(value)
        self.assertIn("P3064/S2014", MD.read_text(encoding="utf-8"))
        self.assertIn("P3064/S2014", (ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md").read_text(encoding="utf-8"))
        self.assertIn("P3064/S2014", (ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md").read_text(encoding="utf-8"))
        self.assertIn("P3064", (REPO / "AGENTS.md").read_text(encoding="utf-8"))

if __name__ == "__main__":
    unittest.main()
