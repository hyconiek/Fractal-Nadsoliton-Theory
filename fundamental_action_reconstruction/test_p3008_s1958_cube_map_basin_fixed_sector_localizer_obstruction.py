import json, subprocess, sys, unittest
from pathlib import Path
ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
SCRIPT = ROOT / "p3008_s1958_cube_map_basin_fixed_sector_localizer_obstruction.py"
OUT = ROOT / "generated" / "p3008_s1958_cube_map_basin_fixed_sector_localizer_obstruction.json"
MD = ROOT / "generated" / "p3008_s1958_cube_map_basin_fixed_sector_localizer_obstruction.md"

class P3008CubeMapBasinFixedSectorLocalizerObstructionTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        subprocess.run([sys.executable, str(SCRIPT)], cwd=ROOT, check=True, stdout=subprocess.PIPE, text=True)
        cls.payload = json.loads(OUT.read_text(encoding="utf-8"))

    def test_status_and_input(self):
        self.assertEqual(self.payload["status"], "P3008_CUBE_MAP_BASIN_FIXED_SECTOR_LOCALIZER_OBSTRUCTION_BOUNDED_NO_GO")
        self.assertIsNotNone(self.payload["input_hashes"]["P3007"])

    def test_localizer_certificate(self):
        cert = self.payload["localizer_certificate"]
        self.assertEqual(cert["fixed_residue_count"], 9)
        self.assertEqual(cert["basin_count"], 9)
        self.assertEqual(cert["residue_row_count"], 12)
        self.assertGreaterEqual(cert["signature_class_count"], 1)
        self.assertEqual(cert["basin_partition_preserving_translations"], [0, 2, 4, 6, 8, 10])
        self.assertEqual(cert["basin_partition_preserving_units"], [1, 5, 7, 11])
        self.assertEqual(cert["accepted_nonpremise_localizers"], [])
        self.assertEqual(cert["acceptance_matrix_rows"], 128)
        self.assertEqual(cert["accepted_rows"], 1)

    def test_obligations_and_nonpromotion(self):
        obj = self.payload["constructed_theoretical_objects"]
        self.assertEqual(obj["object"], "CubeMapBasinFixedSectorLocalizer_ObstructionMatrix")
        obligations = {r["obligation"]: r["satisfied"] for r in obj["proof_obligation_rows"]}
        self.assertTrue(obligations["finite_basin_fixed_sector_receivers"])
        self.assertTrue(obligations["fiber_and_signature_partition"])
        self.assertTrue(obligations["symmetry_nonlocalizer_witness"])
        self.assertFalse(obligations["strict_provenance_available"])
        self.assertFalse(obligations["nonpremise_physical_sector_theorem"])
        self.assertFalse(obligations["unique_origin_or_sector_localizer"])
        self.assertFalse(obligations["accepted_current_nonpremise_localizer"])
        witness = obj["localizer_witness"]
        self.assertEqual(len(witness["basin_rows"]), 9)
        self.assertEqual(len(witness["residue_rows"]), 12)
        for value in self.payload["decision"]["negative_export_flags"].values():
            self.assertFalse(value)

    def test_docs(self):
        self.assertIn("P3008/S1958", MD.read_text(encoding="utf-8"))
        self.assertIn("P3008/S1958", (ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md").read_text(encoding="utf-8"))
        self.assertIn("P3008/S1958", (ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md").read_text(encoding="utf-8"))
        self.assertIn("P3008", (REPO / "AGENTS.md").read_text(encoding="utf-8"))

if __name__ == "__main__":
    unittest.main()
