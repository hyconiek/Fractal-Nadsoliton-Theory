import json, subprocess, sys, unittest
from pathlib import Path
ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
SCRIPT = ROOT / "p3033_s1983_z12_block_coarse_graining_classical_limit_obstruction.py"
OUT = ROOT / "generated" / "p3033_s1983_z12_block_coarse_graining_classical_limit_obstruction.json"
MD = ROOT / "generated" / "p3033_s1983_z12_block_coarse_graining_classical_limit_obstruction.md"

class P3033Z12BlockCoarseGrainingTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        subprocess.run([sys.executable, str(SCRIPT)], cwd=ROOT, check=True, stdout=subprocess.PIPE, text=True)
        cls.payload = json.loads(OUT.read_text(encoding="utf-8"))

    def test_status_and_inputs(self):
        self.assertEqual(self.payload["status"], "P3033_Z12_BLOCK_COARSE_GRAINING_CLASSICAL_LIMIT_OBSTRUCTION_NO_CLASSICAL_EXPORT")
        self.assertIsNotNone(self.payload["input_hashes"]["P3028"])
        self.assertIsNotNone(self.payload["input_hashes"]["P3032"])

    def test_finite_certificate(self):
        cert = self.payload["finite_certificate"]
        self.assertEqual(cert["block_scales"], 6)
        self.assertGreater(cert["composition_rows"], 0)
        self.assertEqual(cert["composition_rows_passed"], cert["composition_rows"])
        self.assertGreaterEqual(cert["translation_sensitive_nontrivial_scales"], 1)
        self.assertFalse(cert["classical_coarse_graining_limit_exported"])

    def test_obligations_and_exports(self):
        obj = self.payload["constructed_theoretical_objects"]
        self.assertEqual(obj["object"], "Z12BlockCoarseGraining_ClassicalLimitObstructionMatrix")
        obligations = {row["obligation"]: row["satisfied"] for row in obj["proof_obligations"]}
        self.assertTrue(obligations["attacks_single_P3028_foundation_atom"])
        self.assertTrue(obligations["explicit_coarse_graining_operator"])
        self.assertTrue(obligations["finite_rg_composition_law"])
        self.assertFalse(obligations["chart_translation_independent_limit"])
        self.assertFalse(obligations["infinite_refinement_or_continuum_parameter"])
        self.assertFalse(obligations["physical_length_unit_for_scale"])
        self.assertFalse(obligations["observer_independent_classical_readout"])
        for value in self.payload["decision"]["negative_export_flags"].values():
            self.assertFalse(value)

    def test_docs(self):
        self.assertIn("P3033/S1983", MD.read_text(encoding="utf-8"))
        self.assertIn("P3033/S1983", (ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md").read_text(encoding="utf-8"))
        self.assertIn("P3033/S1983", (ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md").read_text(encoding="utf-8"))
        self.assertIn("P3033", (REPO / "AGENTS.md").read_text(encoding="utf-8"))

if __name__ == "__main__":
    unittest.main()
