import json, subprocess, sys, unittest
from pathlib import Path
ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
SCRIPT = ROOT / "p3076_s2026_dirichlet_spectral_dispersion_audit.py"
OUT = ROOT / "generated" / "p3076_s2026_dirichlet_spectral_dispersion_audit.json"
MD = ROOT / "generated" / "p3076_s2026_dirichlet_spectral_dispersion_audit.md"

class P3076DirichletSpectralDispersionAuditTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        subprocess.run([sys.executable, str(SCRIPT)], cwd=ROOT, check=True, stdout=subprocess.PIPE, text=True)
        cls.payload = json.loads(OUT.read_text(encoding="utf-8"))

    def test_status_and_input_hash(self):
        self.assertEqual(self.payload["status"], "P3076_INTERNAL_DIRICHLET_DIFFUSIVE_SPECTRAL_BRANCH_WAVE_LIGHTLIKE_OBSTRUCTION")
        self.assertIsNotNone(self.payload["input_hashes"]["P3075"])

    def test_spectral_counts(self):
        cert = self.payload["finite_certificate"]
        self.assertEqual(cert["content_grep_lanes"], 3)
        self.assertGreater(cert["content_grep_hits"], 0)
        self.assertEqual(cert["p3075_local_dirichlet_accepted_rows"], 96)
        self.assertEqual(cert["z12_modes"], 12)
        self.assertEqual(cert["nonconstant_modes"], 11)
        self.assertEqual(cert["fractional_steps"], 3)
        self.assertEqual(cert["amplification_rows"], 36)
        self.assertEqual(cert["nonexpansive_amplification_rows"], 36)
        self.assertEqual(cert["strictly_contracting_nonconstant_amplification_rows"], 33)
        self.assertEqual(cert["oscillatory_phase_rows"], 0)
        self.assertEqual(cert["accepted_lightlike_branch_rows"], 0)
        self.assertEqual(cert["satisfied_proof_obligations"], 4)

    def test_mode_spectrum_labels(self):
        modes = self.payload["constructed_theoretical_objects"]["mode_spectrum_rows"]
        by_j = {row["mode_j"]: row for row in modes}
        self.assertEqual(by_j[0]["laplacian_eigenvalue_exact"], "0")
        self.assertEqual(by_j[1]["laplacian_eigenvalue_exact"], "2-sqrt(3)")
        self.assertEqual(by_j[6]["laplacian_eigenvalue_exact"], "4")
        self.assertEqual(by_j[11]["laplacian_eigenvalue_exact"], "2-sqrt(3)")
        self.assertFalse(any(row["accepted_as_lightlike_branch"] for row in modes))

    def test_flags_recommendation_and_docs(self):
        for value in self.payload["decision"]["negative_export_flags"].values():
            self.assertFalse(value)
        self.assertTrue(self.payload["decision"]["positive_scoped_flags"]["internal_diffusive_smoothing_branch_exported"])
        self.assertIn("second-order lift obstruction table", self.payload["decision"]["next_honest_step"])
        self.assertIn("P3076/S2026", MD.read_text(encoding="utf-8"))
        self.assertIn("P3076/S2026", (ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md").read_text(encoding="utf-8"))
        self.assertIn("P3076/S2026", (ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md").read_text(encoding="utf-8"))
        self.assertIn("P3076", (REPO / "AGENTS.md").read_text(encoding="utf-8"))

if __name__ == "__main__":
    unittest.main()
