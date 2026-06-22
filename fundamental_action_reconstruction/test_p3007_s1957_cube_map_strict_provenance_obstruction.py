import json, subprocess, sys, unittest
from pathlib import Path
ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
SCRIPT = ROOT / "p3007_s1957_cube_map_strict_provenance_obstruction.py"
OUT = ROOT / "generated" / "p3007_s1957_cube_map_strict_provenance_obstruction.json"
MD = ROOT / "generated" / "p3007_s1957_cube_map_strict_provenance_obstruction.md"

class P3007CubeMapStrictProvenanceObstructionTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        subprocess.run([sys.executable, str(SCRIPT)], cwd=ROOT, check=True, stdout=subprocess.PIPE, text=True)
        cls.payload = json.loads(OUT.read_text(encoding="utf-8"))

    def test_status_and_input(self):
        self.assertEqual(self.payload["status"], "P3007_CUBE_MAP_STRICT_PROVENANCE_OBSTRUCTION_BOUNDED_NO_GO")
        self.assertIsNotNone(self.payload["input_hashes"]["P3006"])

    def test_provenance_certificate(self):
        cert = self.payload["provenance_certificate"]
        self.assertEqual(cert["residue_count"], 12)
        self.assertEqual(cert["fixed_residues"], [0, 1, 3, 4, 5, 7, 8, 9, 11])
        self.assertEqual(cert["moved_residues"], [2, 6, 10])
        self.assertEqual(cert["unit_residues"], [1, 5, 7, 11])
        self.assertTrue(cert["all_defects_zero_mod3"])
        self.assertEqual(cert["mod4_defect_carriers"], [2, 6, 10])
        self.assertEqual(cert["strict_provenance_receivers_accepted"], [])
        self.assertEqual(cert["acceptance_matrix_rows"], 512)
        self.assertEqual(cert["accepted_rows"], 1)

    def test_receiver_matrix_and_nonpromotion(self):
        obj = self.payload["constructed_theoretical_objects"]
        self.assertEqual(obj["object"], "CubeMapStrictProvenance_ReceiverObstructionMatrix")
        self.assertEqual(len(obj["crt_defect_rows"]), 12)
        receivers = {r["receiver"]: r["satisfied"] for r in obj["provenance_receiver_matrix"]}
        self.assertTrue(receivers["finite_z12_power_map"])
        self.assertTrue(receivers["crt_mod3_identity"])
        self.assertTrue(receivers["unit_fixed_sector"])
        self.assertTrue(receivers["nonidentity_defect_witness"])
        self.assertTrue(receivers["crt_mod4_even_defect_witness"])
        self.assertFalse(receivers["strict_nadsoliton_cube_law_export"])
        self.assertFalse(receivers["apd_or_phase_provenance_export"])
        self.assertFalse(receivers["damping_compression_provenance_export"])
        self.assertFalse(receivers["accepted_current_strict_provenance"])
        for value in self.payload["decision"]["negative_export_flags"].values():
            self.assertFalse(value)

    def test_docs(self):
        self.assertIn("P3007/S1957", MD.read_text(encoding="utf-8"))
        self.assertIn("P3007/S1957", (ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md").read_text(encoding="utf-8"))
        self.assertIn("P3007/S1957", (ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md").read_text(encoding="utf-8"))
        self.assertIn("P3007", (REPO / "AGENTS.md").read_text(encoding="utf-8"))

if __name__ == "__main__":
    unittest.main()
