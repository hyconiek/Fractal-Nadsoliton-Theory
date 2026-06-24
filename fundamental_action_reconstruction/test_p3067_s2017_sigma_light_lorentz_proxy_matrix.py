import json, subprocess, sys, unittest
from pathlib import Path
ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
SCRIPT = ROOT / "p3067_s2017_sigma_light_lorentz_proxy_matrix.py"
OUT = ROOT / "generated" / "p3067_s2017_sigma_light_lorentz_proxy_matrix.json"
MD = ROOT / "generated" / "p3067_s2017_sigma_light_lorentz_proxy_matrix.md"

class P3067SigmaLightLorentzProxyTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        subprocess.run([sys.executable, str(SCRIPT)], cwd=ROOT, check=True, stdout=subprocess.PIPE, text=True)
        cls.payload = json.loads(OUT.read_text(encoding="utf-8"))

    def test_status_and_input_hash(self):
        self.assertEqual(self.payload["status"], "P3067_SIGMA_LIGHT_LORENTZ_PROXY_MATRIX_CONDITIONAL_NO_OBSERVED_LIGHT_CLOSURE")
        self.assertIsNotNone(self.payload["input_hashes"]["P3066"])

    def test_content_first_and_finite_counts(self):
        cert = self.payload["finite_certificate"]
        self.assertEqual(cert["content_grep_lanes"], 3)
        self.assertGreater(cert["content_grep_hits"], 0)
        self.assertEqual(cert["sigma_branches"], 2)
        self.assertEqual(cert["sampled_boosts"], 4)
        self.assertEqual(cert["lorentz_proxy_rows"], 8)
        self.assertEqual(cert["proxy_null_covariance_pass_rows"], 8)
        self.assertEqual(cert["strict_lorentz_closure_rows"], 0)
        self.assertEqual(cert["p3066_accepted_physics_obligation_rows"], 0)
        self.assertEqual(cert["satisfied_proof_obligations"], 3)

    def test_constructed_transition_proxy_and_rows(self):
        obj = self.payload["constructed_theoretical_objects"]["sigma_conditioned_nadsoliton_to_light_transition_proxy"]
        self.assertEqual(obj["object"], "SigmaConditionedNadsolitonLightTransitionProxy")
        self.assertIn("k_sigma=(1,sigma)", obj["definition"])
        rows = self.payload["constructed_theoretical_objects"]["finite_lorentz_proxy_matrix"]
        self.assertEqual(len(rows), 8)
        self.assertTrue(all(row["proxy_null_covariance_pass"] for row in rows))
        self.assertTrue(all(not row["strict_lorentz_closure_exported"] for row in rows))

    def test_flags_recommendation_and_docs(self):
        for value in self.payload["decision"]["negative_export_flags"].values():
            self.assertFalse(value)
        self.assertTrue(self.payload["decision"]["positive_scoped_flags"]["finite_null_covariance_proxy_passed"])
        self.assertIn("spacetime embedding", self.payload["decision"]["next_honest_step"])
        self.assertIn("P3067/S2017", MD.read_text(encoding="utf-8"))
        self.assertIn("P3067/S2017", (ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md").read_text(encoding="utf-8"))
        self.assertIn("P3067/S2017", (ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md").read_text(encoding="utf-8"))
        self.assertIn("P3067", (REPO / "AGENTS.md").read_text(encoding="utf-8"))

if __name__ == "__main__":
    unittest.main()
