import json, subprocess, sys, unittest
from pathlib import Path
ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
SCRIPT = ROOT / "p3038_s1988_viscous_retarded_chiral_selector_candidate_obstruction.py"
OUT = ROOT / "generated" / "p3038_s1988_viscous_retarded_chiral_selector_candidate_obstruction.json"
MD = ROOT / "generated" / "p3038_s1988_viscous_retarded_chiral_selector_candidate_obstruction.md"

class P3038ViscousRetardedChiralSelectorCandidateTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        subprocess.run([sys.executable, str(SCRIPT)], cwd=ROOT, check=True, stdout=subprocess.PIPE, text=True)
        cls.payload = json.loads(OUT.read_text(encoding="utf-8"))

    def test_status_and_input(self):
        self.assertEqual(self.payload["status"], "P3038_VISCOUS_RETARDED_CHIRAL_SELECTOR_CANDIDATE_OBSTRUCTION_NO_SELECTOR_EXPORT")
        self.assertIsNotNone(self.payload["input_hashes"]["P3037"])

    def test_finite_branch_score(self):
        cert = self.payload["finite_certificate"]
        self.assertEqual(cert["operator_rows"], 4)
        self.assertEqual(cert["finite_nonzero_operator_rows"], 4)
        self.assertEqual(cert["accepted_operator_source_rows"], 0)
        self.assertEqual(cert["p3037_features_present"], 6)
        self.assertTrue(cert["branch_score_nonzero"])
        self.assertGreater(cert["branch_separation_abs"], 0.0)
        self.assertFalse(cert["strict_selector_mechanism_exported"])

    def test_obligations_and_exports(self):
        obj = self.payload["constructed_theoretical_objects"]
        self.assertEqual(obj["object"], "ViscousRetardedChiralReadoutSelector_CandidateObstructionMatrix")
        scores = obj["branch_scores"]
        self.assertAlmostEqual(scores["plus"], -scores["minus"], places=12)
        obligations = {row["obligation"]: row["satisfied"] for row in obj["proof_obligations"]}
        self.assertTrue(obligations["integrated_candidate_operator_constructed"])
        self.assertTrue(obligations["all_p3037_hint_features_present"])
        self.assertTrue(obligations["finite_branch_score_nonzero"])
        self.assertTrue(obligations["branch_separation_computable"])
        self.assertFalse(obligations["nonpremise_chiral_sign_localizer"])
        self.assertFalse(obligations["sourced_retardation_path_anisotropy"])
        self.assertFalse(obligations["physical_unit_readout_coupling"])
        self.assertFalse(obligations["nonproxy_variational_or_readout_theorem"])
        for value in self.payload["decision"]["negative_export_flags"].values():
            self.assertFalse(value)

    def test_docs(self):
        self.assertIn("P3038/S1988", MD.read_text(encoding="utf-8"))
        self.assertIn("P3038/S1988", (ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md").read_text(encoding="utf-8"))
        self.assertIn("P3038/S1988", (ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md").read_text(encoding="utf-8"))
        self.assertIn("P3038", (REPO / "AGENTS.md").read_text(encoding="utf-8"))

if __name__ == "__main__":
    unittest.main()
