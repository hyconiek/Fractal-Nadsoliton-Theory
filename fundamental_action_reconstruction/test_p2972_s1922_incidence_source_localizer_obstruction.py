import json, subprocess, sys, unittest
from pathlib import Path
ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
SCRIPT = ROOT / "p2972_s1922_incidence_source_localizer_obstruction.py"
OUT = ROOT / "generated" / "p2972_s1922_incidence_source_localizer_obstruction.json"
MD = ROOT / "generated" / "p2972_s1922_incidence_source_localizer_obstruction.md"

class P2972IncidenceSourceLocalizerObstructionTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        subprocess.run([sys.executable, str(SCRIPT)], cwd=ROOT, check=True, stdout=subprocess.PIPE, text=True)
        cls.payload = json.loads(OUT.read_text(encoding="utf-8"))

    def test_status_and_input(self):
        self.assertEqual(self.payload["status"], "P2972_INCIDENCE_SOURCE_LOCALIZER_OBSTRUCTION_NO_STRICT_EXPORT")
        self.assertIsNotNone(self.payload["input_hashes"]["P2971"])

    def test_localizer_certificate(self):
        cert = self.payload["localizer_certificate"]
        self.assertEqual(cert["automorphism_count"], 6)
        self.assertEqual(cert["orbit_count"], 3)
        self.assertEqual(cert["invariant_subset_count"], 8)
        self.assertEqual(cert["candidate_count"], 6)
        self.assertEqual(cert["accepted_current_strict_localizers"], [])
        self.assertEqual(cert["acceptance_matrix_rows"], 64)
        self.assertEqual(cert["accepted_rows"], 1)

    def test_obligations_and_rows(self):
        obligations = {r["obligation"]: r["satisfied"] for r in self.payload["constructed_theoretical_objects"]["proof_obligation_rows"]}
        self.assertTrue(obligations["slot_orbits_computed"])
        self.assertTrue(obligations["invariant_localizer_candidate_exists"])
        self.assertTrue(obligations["whole_object_candidate_exists"])
        self.assertFalse(obligations["strict_source_theorem_exported"])
        self.assertFalse(obligations["unit_coupling_exported"])
        self.assertFalse(obligations["accepted_current_strict_localizer"])
        rows = {r["candidate"]: r for r in self.payload["constructed_theoretical_objects"]["localizer_rows"]}
        self.assertTrue(rows["whole_incidence_identity_localizer"]["selects_whole_incidence_object"])
        self.assertFalse(rows["aggregate_coordinate_order_localizer"]["automorphism_invariant"])
        self.assertFalse(rows["completed_strict_incidence_source_localizer_schema"]["accepted_current_strict_localizer"])

    def test_docs_and_nonpromotion(self):
        for value in self.payload["decision"]["negative_export_flags"].values():
            self.assertFalse(value)
        self.assertIn("P2972/S1922", MD.read_text(encoding="utf-8"))
        self.assertIn("P2972/S1922", (ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md").read_text(encoding="utf-8"))
        self.assertIn("P2972/S1922", (ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md").read_text(encoding="utf-8"))
        self.assertIn("P2972", (REPO / "AGENTS.md").read_text(encoding="utf-8"))

if __name__ == "__main__":
    unittest.main()
