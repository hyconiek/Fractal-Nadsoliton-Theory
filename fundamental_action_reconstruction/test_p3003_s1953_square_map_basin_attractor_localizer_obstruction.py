import json, subprocess, sys, unittest
from pathlib import Path
ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
SCRIPT = ROOT / "p3003_s1953_square_map_basin_attractor_localizer_obstruction.py"
OUT = ROOT / "generated" / "p3003_s1953_square_map_basin_attractor_localizer_obstruction.json"
MD = ROOT / "generated" / "p3003_s1953_square_map_basin_attractor_localizer_obstruction.md"

class P3003SquareMapBasinAttractorLocalizerObstructionTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        subprocess.run([sys.executable, str(SCRIPT)], cwd=ROOT, check=True, stdout=subprocess.PIPE, text=True)
        cls.payload = json.loads(OUT.read_text(encoding="utf-8"))

    def test_status_and_input(self):
        self.assertEqual(self.payload["status"], "P3003_SQUARE_MAP_BASIN_ATTRACTOR_LOCALIZER_OBSTRUCTION_BOUNDED_NO_GO")
        self.assertIsNotNone(self.payload["input_hashes"]["P3002"])

    def test_localizer_certificate(self):
        cert = self.payload["localizer_certificate"]
        self.assertEqual(cert["basin_row_count"], 4)
        self.assertEqual(cert["residue_row_count"], 12)
        self.assertGreaterEqual(cert["signature_class_count"], 4)
        self.assertIn(0, cert["singleton_signature_classes"])
        self.assertEqual(cert["basin_partition_preserving_translations"], [0, 3, 6, 9])
        self.assertEqual(cert["accepted_nonpremise_localizer_count"], 0)
        self.assertEqual(cert["acceptance_matrix_rows"], 64)
        self.assertEqual(cert["accepted_rows"], 1)

    def test_obligations_and_nonpromotion(self):
        obj = self.payload["constructed_theoretical_objects"]
        obligations = {r["obligation"]: r["satisfied"] for r in obj["proof_obligation_rows"]}
        self.assertTrue(obligations["finite_basin_attractor_receivers"])
        self.assertTrue(obligations["signature_partition_computed"])
        self.assertTrue(obligations["translation_nonlocalizer_witness"])
        self.assertFalse(obligations["strict_provenance_available"])
        self.assertFalse(obligations["nonpremise_physical_sector_theorem"])
        self.assertFalse(obligations["accepted_current_nonpremise_localizer"])
        witness = obj["localizer_witness"]
        self.assertEqual(len(witness["basin_rows"]), 4)
        self.assertEqual(len(witness["residue_rows"]), 12)
        for value in self.payload["decision"]["negative_export_flags"].values():
            self.assertFalse(value)

    def test_docs(self):
        self.assertIn("P3003/S1953", MD.read_text(encoding="utf-8"))
        self.assertIn("P3003/S1953", (ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md").read_text(encoding="utf-8"))
        self.assertIn("P3003/S1953", (ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md").read_text(encoding="utf-8"))
        self.assertIn("P3003", (REPO / "AGENTS.md").read_text(encoding="utf-8"))

if __name__ == "__main__":
    unittest.main()
