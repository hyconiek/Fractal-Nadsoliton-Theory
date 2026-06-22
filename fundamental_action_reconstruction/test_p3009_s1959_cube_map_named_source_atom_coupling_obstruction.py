import json, subprocess, sys, unittest
from pathlib import Path
ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
SCRIPT = ROOT / "p3009_s1959_cube_map_named_source_atom_coupling_obstruction.py"
OUT = ROOT / "generated" / "p3009_s1959_cube_map_named_source_atom_coupling_obstruction.json"
MD = ROOT / "generated" / "p3009_s1959_cube_map_named_source_atom_coupling_obstruction.md"

class P3009CubeMapNamedSourceAtomCouplingObstructionTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        subprocess.run([sys.executable, str(SCRIPT)], cwd=ROOT, check=True, stdout=subprocess.PIPE, text=True)
        cls.payload = json.loads(OUT.read_text(encoding="utf-8"))

    def test_status_and_input(self):
        self.assertEqual(self.payload["status"], "P3009_CUBE_MAP_NAMED_SOURCE_ATOM_COUPLING_OBSTRUCTION_BOUNDED_NO_GO")
        self.assertIsNotNone(self.payload["input_hashes"]["P3008"])

    def test_coupling_certificate(self):
        cert = self.payload["coupling_certificate"]
        self.assertEqual(cert["basin_fixed_sector_receiver_count"], 9)
        self.assertEqual(cert["residue_receiver_count"], 12)
        self.assertEqual(cert["total_receiver_count"], 21)
        self.assertEqual(cert["named_source_atom_count"], 4)
        self.assertEqual(cert["coupling_test_count"], 84)
        self.assertEqual(cert["accepted_source_coupling_count"], 0)
        self.assertEqual(cert["acceptance_matrix_rows"], 128)
        self.assertEqual(cert["accepted_rows"], 1)

    def test_obligations_and_nonpromotion(self):
        obj = self.payload["constructed_theoretical_objects"]
        self.assertEqual(obj["object"], "CubeMapNamedSourceAtomCoupling_ObstructionMatrix")
        obligations = {r["obligation"]: r["satisfied"] for r in obj["proof_obligation_rows"]}
        self.assertTrue(obligations["finite_cube_map_atom_matrix"])
        self.assertTrue(obligations["exact_basin_fixed_residue_receivers"])
        self.assertFalse(obligations["strict_cube_map_provenance"])
        self.assertFalse(obligations["accepted_nonpremise_localizer"])
        self.assertFalse(obligations["atom_specific_coupling_theorem"])
        self.assertFalse(obligations["unit_bearing_coupling_coefficient"])
        self.assertFalse(obligations["nonproxy_export"])
        self.assertFalse(obligations["accepted_current_source_coupling"])
        self.assertEqual(len(obj["coupling_witness"]["coupling_rows"]), 84)
        for value in self.payload["decision"]["negative_export_flags"].values():
            self.assertFalse(value)

    def test_docs(self):
        self.assertIn("P3009/S1959", MD.read_text(encoding="utf-8"))
        self.assertIn("P3009/S1959", (ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md").read_text(encoding="utf-8"))
        self.assertIn("P3009/S1959", (ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md").read_text(encoding="utf-8"))
        self.assertIn("P3009", (REPO / "AGENTS.md").read_text(encoding="utf-8"))

if __name__ == "__main__":
    unittest.main()
