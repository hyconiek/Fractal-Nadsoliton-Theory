import json, subprocess, sys, unittest
from pathlib import Path
ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
SCRIPT = ROOT / "p2989_s1939_annihilator_lattice_named_source_atom_coupling_obstruction.py"
OUT = ROOT / "generated" / "p2989_s1939_annihilator_lattice_named_source_atom_coupling_obstruction.json"
MD = ROOT / "generated" / "p2989_s1939_annihilator_lattice_named_source_atom_coupling_obstruction.md"

class P2989AnnihilatorLatticeNamedSourceAtomCouplingObstructionTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        subprocess.run([sys.executable, str(SCRIPT)], cwd=ROOT, check=True, stdout=subprocess.PIPE, text=True)
        cls.payload = json.loads(OUT.read_text(encoding="utf-8"))

    def test_status_and_input(self):
        self.assertEqual(self.payload["status"], "P2989_ANNIHILATOR_LATTICE_NAMED_SOURCE_ATOM_COUPLING_OBSTRUCTION_BOUNDED_NO_GO")
        self.assertIsNotNone(self.payload["input_hashes"]["P2988"])

    def test_coupling_certificate(self):
        cert = self.payload["coupling_certificate"]
        self.assertEqual(cert["row_count"], 6)
        self.assertEqual(cert["coupling_row_count"], 24)
        self.assertEqual(len(cert["named_atoms"]), 4)
        self.assertEqual(cert["accepted_named_atom_couplings"], [])
        self.assertEqual(cert["acceptance_matrix_rows"], 256)
        self.assertEqual(cert["accepted_rows"], 1)

    def test_obligations_and_couplings(self):
        obj = self.payload["constructed_theoretical_objects"]
        obligations = {r["obligation"]: r["satisfied"] for r in obj["proof_obligation_rows"]}
        self.assertTrue(obligations["finite_annihilator_rows_present"])
        self.assertTrue(obligations["named_atom_matrix_constructed"])
        self.assertTrue(obligations["formal_algebraic_receivers_available"])
        self.assertFalse(obligations["strict_provenance_available"])
        self.assertFalse(obligations["nonpremise_source_localizer_available"])
        self.assertFalse(obligations["atom_specific_source_theorem"])
        self.assertFalse(obligations["accepted_current_named_atom_coupling"])
        rows = obj["coupling_witness"]["coupling_rows"]
        self.assertEqual(sum(len(r["atom_couplings"]) for r in rows), 24)
        self.assertFalse(any(c["accepted_named_atom_coupling"] for r in rows for c in r["atom_couplings"]))

    def test_docs_and_nonpromotion(self):
        for value in self.payload["decision"]["negative_export_flags"].values():
            self.assertFalse(value)
        self.assertIn("P2989/S1939", MD.read_text(encoding="utf-8"))
        self.assertIn("P2989/S1939", (ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md").read_text(encoding="utf-8"))
        self.assertIn("P2989/S1939", (ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md").read_text(encoding="utf-8"))
        self.assertIn("P2989", (REPO / "AGENTS.md").read_text(encoding="utf-8"))

if __name__ == "__main__":
    unittest.main()
