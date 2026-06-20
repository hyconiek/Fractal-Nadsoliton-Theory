import json, subprocess, sys, unittest
from pathlib import Path
ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
SCRIPT = ROOT / "p2986_s1936_z12_annihilator_lattice_source_candidate_obstruction.py"
OUT = ROOT / "generated" / "p2986_s1936_z12_annihilator_lattice_source_candidate_obstruction.json"
MD = ROOT / "generated" / "p2986_s1936_z12_annihilator_lattice_source_candidate_obstruction.md"

class P2986Z12AnnihilatorLatticeSourceCandidateObstructionTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        subprocess.run([sys.executable, str(SCRIPT)], cwd=ROOT, check=True, stdout=subprocess.PIPE, text=True)
        cls.payload = json.loads(OUT.read_text(encoding="utf-8"))

    def test_status_and_input(self):
        self.assertEqual(self.payload["status"], "P2986_Z12_ANNIHILATOR_LATTICE_SOURCE_CANDIDATE_OBSTRUCTION_BOUNDED_NO_GO")
        self.assertIsNotNone(self.payload["input_hashes"]["P2985"])

    def test_annihilator_certificate(self):
        cert = self.payload["annihilator_certificate"]
        self.assertEqual(cert["modulus"], 12)
        self.assertEqual(cert["ideal_count"], 6)
        self.assertTrue(cert["all_annihilator_products_zero"])
        self.assertTrue(cert["all_double_annihilator_exact"])
        self.assertTrue(cert["all_order_reversing_checks_pass"])
        self.assertEqual(cert["accepted_current_source_candidates"], [])
        self.assertEqual(cert["acceptance_matrix_rows"], 256)
        self.assertEqual(cert["accepted_rows"], 1)

    def test_obligations_and_rows(self):
        obj = self.payload["constructed_theoretical_objects"]
        obligations = {r["obligation"]: r["satisfied"] for r in obj["proof_obligation_rows"]}
        self.assertTrue(obligations["finite_ideal_lattice_enumerated"])
        self.assertTrue(obligations["annihilator_products_zero"])
        self.assertTrue(obligations["double_annihilator_closure"])
        self.assertTrue(obligations["order_reversing_lattice_witness"])
        self.assertFalse(obligations["strict_nadsoliton_provenance"])
        self.assertFalse(obligations["nonpremise_source_localizer"])
        self.assertFalse(obligations["named_source_atom_coupling"])
        self.assertFalse(obligations["positive_unit_or_measure_source"])
        self.assertFalse(obligations["unit_bearing_action_installation"])
        rows = obj["annihilator_lattice_witness"]["annihilator_rows"]
        self.assertEqual(len(rows), 6)
        self.assertTrue(any(r["ideal"] == [0] and r["annihilator"] == list(range(12)) for r in rows))
        self.assertTrue(any(r["ideal"] == list(range(12)) and r["annihilator"] == [0] for r in rows))

    def test_docs_and_nonpromotion(self):
        for value in self.payload["decision"]["negative_export_flags"].values():
            self.assertFalse(value)
        self.assertIn("P2986/S1936", MD.read_text(encoding="utf-8"))
        self.assertIn("P2986/S1936", (ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md").read_text(encoding="utf-8"))
        self.assertIn("P2986/S1936", (ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md").read_text(encoding="utf-8"))
        self.assertIn("P2986", (REPO / "AGENTS.md").read_text(encoding="utf-8"))

if __name__ == "__main__":
    unittest.main()
