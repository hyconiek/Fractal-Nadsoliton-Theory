import json, subprocess, sys, unittest
from pathlib import Path
ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
SCRIPT = ROOT / "p3040_s1990_retarded_path_anisotropy_source_obstruction.py"
OUT = ROOT / "generated" / "p3040_s1990_retarded_path_anisotropy_source_obstruction.json"
MD = ROOT / "generated" / "p3040_s1990_retarded_path_anisotropy_source_obstruction.md"

class P3040RetardedPathAnisotropySourceTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        subprocess.run([sys.executable, str(SCRIPT)], cwd=ROOT, check=True, stdout=subprocess.PIPE, text=True)
        cls.payload = json.loads(OUT.read_text(encoding="utf-8"))

    def test_status_and_inputs(self):
        self.assertEqual(self.payload["status"], "P3040_RETARDED_PATH_ANISOTROPY_SOURCE_OBSTRUCTION_NO_SOURCE_EXPORT")
        self.assertIsNotNone(self.payload["input_hashes"]["P3038"])
        self.assertIsNotNone(self.payload["input_hashes"]["P3039"])

    def test_finite_certificate(self):
        cert = self.payload["finite_certificate"]
        self.assertEqual(cert["receiver_rows"], 4)
        self.assertEqual(cert["finite_nonzero_rows"], 4)
        self.assertEqual(cert["accepted_path_anisotropy_source_rows"], 0)
        self.assertEqual(cert["rho_grid_rows"], 12)
        self.assertTrue(cert["path_sign_flip_verified"])
        self.assertEqual(cert["satisfied_proof_obligations"], 4)
        self.assertFalse(cert["sourced_retardation_path_anisotropy_exported"])

    def test_obligations_and_exports(self):
        obj = self.payload["constructed_theoretical_objects"]
        self.assertEqual(obj["object"], "RetardedPathAnisotropySource_ObstructionMatrix")
        obligations = {row["obligation"]: row["satisfied"] for row in obj["proof_obligations"]}
        self.assertTrue(obligations["exact_p3038_path_anisotropy_premise_targeted"])
        self.assertTrue(obligations["finite_retarded_split_torsor_constructed"])
        self.assertTrue(obligations["nonzero_candidate_split_verified"])
        self.assertTrue(obligations["path_sign_flip_verified"])
        self.assertFalse(obligations["nonpremise_positive_rho_source"])
        self.assertFalse(obligations["strict_parallel_perpendicular_path_geometry"])
        self.assertFalse(obligations["chart_independent_path_localizer"])
        self.assertFalse(obligations["coupling_back_to_p3038_as_source"])
        for value in self.payload["decision"]["negative_export_flags"].values():
            self.assertFalse(value)

    def test_docs(self):
        self.assertIn("P3040/S1990", MD.read_text(encoding="utf-8"))
        self.assertIn("P3040/S1990", (ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md").read_text(encoding="utf-8"))
        self.assertIn("P3040/S1990", (ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md").read_text(encoding="utf-8"))
        self.assertIn("P3040", (REPO / "AGENTS.md").read_text(encoding="utf-8"))

if __name__ == "__main__":
    unittest.main()
