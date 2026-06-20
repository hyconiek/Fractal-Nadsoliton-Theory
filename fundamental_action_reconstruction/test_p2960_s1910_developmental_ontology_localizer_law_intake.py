import json, subprocess, sys, unittest
from pathlib import Path
ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
SCRIPT = ROOT / "p2960_s1910_developmental_ontology_localizer_law_intake.py"
OUT = ROOT / "generated" / "p2960_s1910_developmental_ontology_localizer_law_intake.json"
MD = ROOT / "generated" / "p2960_s1910_developmental_ontology_localizer_law_intake.md"

class P2960DevelopmentalOntologyLocalizerLawIntakeTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        subprocess.run([sys.executable, str(SCRIPT)], cwd=ROOT, check=True, stdout=subprocess.PIPE, text=True)
        cls.payload = json.loads(OUT.read_text(encoding="utf-8"))

    def test_status_and_input(self):
        self.assertEqual(self.payload["status"], "P2960_DEVELOPMENTAL_ONTOLOGY_LOCALIZER_LAW_INTAKE_NO_STRICT_EXPORT")
        self.assertIsNotNone(self.payload["input_hashes"]["P2959"])

    def test_developmental_candidates(self):
        cert = self.payload["intake_certificate"]
        self.assertEqual(cert["candidate_law_count"], 4)
        self.assertIn("developmental_fractal_self_balance_minimal_positive_source", cert["new_developmental_candidates_selecting_target"])
        self.assertIn("developmental_minimax_source_amplitude_quotient", cert["new_developmental_candidates_selecting_target"])
        self.assertTrue(cert["current_artifact_accepts_as_developmental_work_item"])
        self.assertFalse(cert["current_artifact_accepts_as_strict_closure"])

    def test_acceptance_matrix_and_obligations(self):
        objs = self.payload["constructed_theoretical_objects"]
        self.assertEqual(len(objs["finite_acceptance_matrix"]), 64)
        self.assertEqual(sum(1 for r in objs["finite_acceptance_matrix"] if r["accepts_as_strict_closure"]), 1)
        self.assertEqual(sum(1 for r in objs["finite_acceptance_matrix"] if r["accepts_as_developmental_work_item"]), 8)
        obligations = {r["obligation"]: r["satisfied"] for r in objs["acceptance_obligation_rows"]}
        self.assertTrue(obligations["content_replay_excluded"])
        self.assertTrue(obligations["developmental_theory_candidate_constructed"])
        self.assertTrue(obligations["ontology_respects_nadsoliton_primordiality"])
        self.assertFalse(obligations["nonconventional_source_derivation_exported"])
        self.assertFalse(obligations["canonical_scale_quotient_exported"])
        self.assertFalse(obligations["unit_bearing_nonproxy_coupling_exported"])

    def test_docs_and_nonpromotion(self):
        for value in self.payload["decision"]["negative_export_flags"].values():
            self.assertFalse(value)
        self.assertIn("P2960/S1910", MD.read_text(encoding="utf-8"))
        self.assertIn("P2960/S1910", (ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md").read_text(encoding="utf-8"))
        self.assertIn("P2960/S1910", (ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md").read_text(encoding="utf-8"))
        self.assertIn("P2960", (REPO / "AGENTS.md").read_text(encoding="utf-8"))

if __name__ == "__main__":
    unittest.main()
