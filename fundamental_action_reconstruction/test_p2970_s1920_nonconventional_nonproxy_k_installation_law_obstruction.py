import json, subprocess, sys, unittest
from pathlib import Path
ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
SCRIPT = ROOT / "p2970_s1920_nonconventional_nonproxy_k_installation_law_obstruction.py"
OUT = ROOT / "generated" / "p2970_s1920_nonconventional_nonproxy_k_installation_law_obstruction.json"
MD = ROOT / "generated" / "p2970_s1920_nonconventional_nonproxy_k_installation_law_obstruction.md"

class P2970NonconventionalNonproxyKInstallationLawObstructionTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        subprocess.run([sys.executable, str(SCRIPT)], cwd=ROOT, check=True, stdout=subprocess.PIPE, text=True)
        cls.payload = json.loads(OUT.read_text(encoding="utf-8"))

    def test_status_and_inputs(self):
        self.assertEqual(self.payload["status"], "P2970_NONCONVENTIONAL_NONPROXY_K_INSTALLATION_LAW_OBSTRUCTION_NO_STRICT_EXPORT")
        self.assertIsNotNone(self.payload["input_hashes"]["P2968"])
        self.assertIsNotNone(self.payload["input_hashes"]["P2969"])

    def test_installation_certificate(self):
        cert = self.payload["installation_certificate"]
        self.assertEqual(cert["distinct_exponent_count"], 24)
        self.assertEqual(cert["candidate_count"], 7)
        self.assertEqual(cert["accepted_current_strict_installations"], [])
        self.assertEqual(cert["acceptance_matrix_rows"], 128)
        self.assertEqual(cert["accepted_rows"], 1)

    def test_obligations_and_rows(self):
        obligations = {r["obligation"]: r["satisfied"] for r in self.payload["constructed_theoretical_objects"]["proof_obligation_rows"]}
        self.assertTrue(obligations["P2968_exponent_set_reused"])
        self.assertTrue(obligations["unique_k_candidate_exists"])
        self.assertTrue(obligations["nonconventional_current_candidate_exists"])
        self.assertFalse(obligations["nonproxy_action_density_coupled"])
        self.assertTrue(obligations["no_blocked_lane_replay"])
        self.assertFalse(obligations["accepted_current_strict_k_installation"])
        rows = {r["candidate"]: r for r in self.payload["constructed_theoretical_objects"]["installation_law_rows"]}
        self.assertEqual(rows["minimal_absolute_exponent_installation"]["selected_k"], ["0/1"])
        self.assertFalse(rows["nonproxy_euler_stationarity_placeholder"]["unique_k"])
        self.assertFalse(rows["completed_strict_k_installation_law_schema"]["current_artifact_available"])

    def test_docs_and_nonpromotion(self):
        for value in self.payload["decision"]["negative_export_flags"].values():
            self.assertFalse(value)
        self.assertIn("P2970/S1920", MD.read_text(encoding="utf-8"))
        self.assertIn("P2970/S1920", (ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md").read_text(encoding="utf-8"))
        self.assertIn("P2970/S1920", (ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md").read_text(encoding="utf-8"))
        self.assertIn("P2970", (REPO / "AGENTS.md").read_text(encoding="utf-8"))

if __name__ == "__main__":
    unittest.main()
