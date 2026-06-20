import json, subprocess, sys, unittest
from pathlib import Path
ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
SCRIPT = ROOT / "p2976_s1926_z12_nilradical_filtration_source_candidate.py"
OUT = ROOT / "generated" / "p2976_s1926_z12_nilradical_filtration_source_candidate.json"
MD = ROOT / "generated" / "p2976_s1926_z12_nilradical_filtration_source_candidate.md"

class P2976Z12NilradicalFiltrationSourceCandidateTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        subprocess.run([sys.executable, str(SCRIPT)], cwd=ROOT, check=True, stdout=subprocess.PIPE, text=True)
        cls.payload = json.loads(OUT.read_text(encoding="utf-8"))

    def test_status_and_input(self):
        self.assertEqual(self.payload["status"], "P2976_Z12_NILRADICAL_FILTRATION_SOURCE_CANDIDATE_DEVELOPMENTAL_NO_STRICT_EXPORT")
        self.assertIsNotNone(self.payload["input_hashes"]["P2975"])

    def test_nilradical_certificate(self):
        cert = self.payload["nilradical_certificate"]
        self.assertEqual(cert["modulus"], 12)
        self.assertEqual(cert["unit_count"], 4)
        self.assertEqual(cert["nilradical"], [0, 6])
        self.assertEqual(cert["nilradical_size"], 2)
        self.assertEqual(cert["nonzero_nilpotents"], [6])
        self.assertEqual(cert["max_first_zero_power"], 2)
        self.assertTrue(cert["all_units_fix_nonzero_nilpotent"])
        self.assertEqual(cert["accepted_current_strict_source_objects"], [])
        self.assertEqual(cert["acceptance_matrix_rows"], 128)
        self.assertEqual(cert["accepted_rows"], 1)

    def test_obligations_and_rows(self):
        obligations = {r["obligation"]: r["satisfied"] for r in self.payload["constructed_theoretical_objects"]["proof_obligation_rows"]}
        self.assertTrue(obligations["new_outside_incidence_lane"])
        self.assertTrue(obligations["not_closed_replay_lane"])
        self.assertTrue(obligations["finite_witness_computed"])
        self.assertFalse(obligations["strict_nadsoliton_provenance_exported"])
        self.assertFalse(obligations["couples_to_named_missing_object"])
        self.assertFalse(obligations["orientation_or_selector_source_exported"])
        self.assertFalse(obligations["action_density_or_variational_installation"])
        rows = {r["candidate"]: r for r in self.payload["constructed_theoretical_objects"]["candidate_rows"]}
        self.assertTrue(rows["Z12_nilradical_filtration_object"]["genuinely_new_outside_incidence_lane"])
        self.assertFalse(rows["completed_strict_nilradical_source_schema"]["accepted_current_strict_source_object"])

    def test_docs_and_nonpromotion(self):
        for value in self.payload["decision"]["negative_export_flags"].values():
            self.assertFalse(value)
        self.assertIn("P2976/S1926", MD.read_text(encoding="utf-8"))
        self.assertIn("P2976/S1926", (ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md").read_text(encoding="utf-8"))
        self.assertIn("P2976/S1926", (ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md").read_text(encoding="utf-8"))
        self.assertIn("P2976", (REPO / "AGENTS.md").read_text(encoding="utf-8"))

if __name__ == "__main__":
    unittest.main()
