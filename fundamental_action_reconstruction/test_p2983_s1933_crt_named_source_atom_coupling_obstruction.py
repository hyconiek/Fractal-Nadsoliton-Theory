import json, subprocess, sys, unittest
from pathlib import Path
ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
SCRIPT = ROOT / "p2983_s1933_crt_named_source_atom_coupling_obstruction.py"
OUT = ROOT / "generated" / "p2983_s1933_crt_named_source_atom_coupling_obstruction.json"
MD = ROOT / "generated" / "p2983_s1933_crt_named_source_atom_coupling_obstruction.md"

class P2983CRTNamedSourceAtomCouplingObstructionTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        subprocess.run([sys.executable, str(SCRIPT)], cwd=ROOT, check=True, stdout=subprocess.PIPE, text=True)
        cls.payload = json.loads(OUT.read_text(encoding="utf-8"))

    def test_status_and_input(self):
        self.assertEqual(self.payload["status"], "P2983_CRT_NAMED_SOURCE_ATOM_COUPLING_OBSTRUCTION_BOUNDED_NO_GO")
        self.assertIsNotNone(self.payload["input_hashes"]["P2982"])

    def test_coupling_certificate(self):
        cert = self.payload["coupling_certificate"]
        self.assertEqual(cert["named_source_atom_count"], 4)
        self.assertEqual(cert["candidate_row_count"], 5)
        self.assertEqual(cert["projectors"], [4, 9])
        self.assertEqual(cert["orthogonal_completion_pair"]["product_mod_12"], 0)
        self.assertEqual(cert["orthogonal_completion_pair"]["sum_mod_12"], 1)
        self.assertEqual(cert["orientation_score_gap"], {"4": 0, "9": 0})
        self.assertEqual(cert["accepted_current_couplings"], [])
        self.assertEqual(cert["acceptance_matrix_rows"], 256)
        self.assertEqual(cert["accepted_rows"], 1)

    def test_obligations_and_rows(self):
        obligations = {r["obligation"]: r["satisfied"] for r in self.payload["constructed_theoretical_objects"]["proof_obligation_rows"]}
        self.assertTrue(obligations["finite_CRT_projector_witness"])
        self.assertTrue(obligations["exactly_named_source_atom_targeted"])
        self.assertFalse(obligations["strict_provenance_available"])
        self.assertFalse(obligations["nonpremise_factor_semantics_available"])
        self.assertFalse(obligations["orientation_sensitive_coupling"])
        self.assertFalse(obligations["positive_scale_or_measure_source"])
        self.assertFalse(obligations["bridge_completion_law"])
        self.assertFalse(obligations["action_variational_installation"])
        self.assertFalse(obligations["accepted_current_coupling"])
        rows = {r["atom"]: r for r in self.payload["constructed_theoretical_objects"]["named_source_atom_rows"]}
        self.assertTrue(rows["selector_orientation_sign_atom"]["finite_crt_witness"])
        self.assertFalse(rows["completed_strict_CRT_named_atom_coupling_schema"]["accepted_current_coupling"])

    def test_docs_and_nonpromotion(self):
        for value in self.payload["decision"]["negative_export_flags"].values():
            self.assertFalse(value)
        self.assertIn("P2983/S1933", MD.read_text(encoding="utf-8"))
        self.assertIn("P2983/S1933", (ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md").read_text(encoding="utf-8"))
        self.assertIn("P2983/S1933", (ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md").read_text(encoding="utf-8"))
        self.assertIn("P2983", (REPO / "AGENTS.md").read_text(encoding="utf-8"))

if __name__ == "__main__":
    unittest.main()
