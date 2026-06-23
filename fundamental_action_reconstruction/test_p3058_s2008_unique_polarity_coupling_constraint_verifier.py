import json, subprocess, sys, unittest
from pathlib import Path
ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
SCRIPT = ROOT / "p3058_s2008_unique_polarity_coupling_constraint_verifier.py"
OUT = ROOT / "generated" / "p3058_s2008_unique_polarity_coupling_constraint_verifier.json"
MD = ROOT / "generated" / "p3058_s2008_unique_polarity_coupling_constraint_verifier.md"

class P3058UniquePolarityCouplingConstraintVerifierTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        subprocess.run([sys.executable, str(SCRIPT)], cwd=ROOT, check=True, stdout=subprocess.PIPE, text=True)
        cls.payload = json.loads(OUT.read_text(encoding="utf-8"))

    def test_status_and_input(self):
        self.assertEqual(self.payload["status"], "P3058_UNIQUE_POLARITY_COUPLING_CONSTRAINT_VERIFIER_BOUNDED_NO_EXPORT")
        self.assertIsNotNone(self.payload["input_hashes"]["P3057"])

    def test_finite_certificate(self):
        cert = self.payload["finite_certificate"]
        self.assertEqual(cert["content_grep_lanes"], 3)
        self.assertGreater(cert["content_grep_hits"], 0)
        self.assertEqual(cert["current_constraints"], 5)
        self.assertEqual(cert["constraint_intersections"], 31)
        self.assertEqual(cert["unique_polarity_intersections"], 0)
        self.assertEqual(cert["all_current_constraints_accepted_count"], 2)
        self.assertEqual(cert["all_current_constraints_accepted_polarities"], ["minus_polarity", "plus_polarity"])
        self.assertEqual(cert["polarity_odd_current_constraints"], 0)
        self.assertEqual(cert["missing_source_atoms_named"], 1)
        self.assertEqual(cert["satisfied_proof_obligations"], 3)

    def test_constructed_object_and_negative_exports(self):
        obj = self.payload["constructed_theoretical_objects"]["unique_polarity_atom_normal_form"]
        self.assertEqual(obj["object"], "UniquePolarityCouplingAtomNormalForm")
        self.assertEqual(obj["target_atom"], "new_unique_polarity_coupling")
        self.assertEqual(obj["carrier_symbol"], "G_selector")
        self.assertEqual(obj["missing_source_atom"]["name"], "strict_polarity_odd_source_law_boundary_condition")
        for row in self.payload["constructed_theoretical_objects"]["constraint_intersections"]:
            self.assertFalse(row["unique"])
            self.assertEqual(row["accepted_count"], 2)
        for value in self.payload["decision"]["negative_export_flags"].values():
            self.assertFalse(value)

    def test_docs(self):
        self.assertIn("P3058/S2008", MD.read_text(encoding="utf-8"))
        self.assertIn("P3058/S2008", (ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md").read_text(encoding="utf-8"))
        self.assertIn("P3058/S2008", (ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md").read_text(encoding="utf-8"))
        self.assertIn("P3058", (REPO / "AGENTS.md").read_text(encoding="utf-8"))

if __name__ == "__main__":
    unittest.main()
