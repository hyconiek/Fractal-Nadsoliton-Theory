import json
import subprocess
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
SCRIPT = ROOT / "p3147_s2097_axiom_lie_smgr_fit_readiness_matrix.py"
OUT = ROOT / "generated" / "p3147_s2097_axiom_lie_smgr_fit_readiness_matrix.json"
MD = ROOT / "generated" / "p3147_s2097_axiom_lie_smgr_fit_readiness_matrix.md"


class P3147AxiomLieSmgrFitReadinessMatrixTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        subprocess.run([sys.executable, str(SCRIPT)], cwd=ROOT, check=True, stdout=subprocess.PIPE, text=True)
        cls.payload = json.loads(OUT.read_text(encoding="utf-8"))

    def test_status_and_finite_counts(self):
        self.assertEqual(self.payload["status"], "P3147_AXIOM_LIE_SMGR_FIT_READINESS_STRONG_CONDITIONAL_NO_STRICT_CLOSURE")
        counts = self.payload["finite_theorem"]["finite_counts"]
        self.assertEqual(counts["rows_audited"], 4)
        self.assertEqual(counts["receiver_or_formal_carrier_rows"], 4)
        self.assertEqual(counts["exact_algebra_or_local_witness_rows"], 4)
        self.assertEqual(counts["conditional_unit_action_axiom_rows"], 4)
        self.assertEqual(counts["full_field_registry_or_metric_bundle_rows"], 0)
        self.assertEqual(counts["nonimported_source_law_rows"], 0)
        self.assertEqual(counts["global_eom_or_bv_brst_closure_rows"], 0)
        self.assertEqual(counts["closed_rows"], 0)
        self.assertGreater(counts["p1960_jacobi_components_checked"], 4000)

    def test_lie_algebra_verdict_is_positive_but_bounded(self):
        rows = {r["row"]: r for r in self.payload["readiness_rows"]}
        self.assertIn("lie_algebra_quality", rows)
        self.assertEqual(rows["lie_algebra_quality"]["readiness_status"], "strong_axiom_branch_potential")
        self.assertIn("good as local algebra", rows["lie_algebra_quality"]["verdict"])
        self.assertIn("global BV/BRST", " ".join(rows["lie_algebra_quality"]["missing_for_strict_closure"]))

    def test_decision_answers_user_questions_without_closure(self):
        decision = self.payload["decision"]
        self.assertIn("Yes", decision["potential_under_axioms"])
        self.assertIn("Good locally/algebraically", decision["lie_algebra_verdict"])
        self.assertIn("conditional receiver/axiom architecture", decision["sm_gr_fit_verdict"])
        self.assertIn("P3148", decision["next_honest_step"])
        self.assertTrue(all(value is False for value in decision["negative_export_flags"].values()))

    def test_docs_updated(self):
        self.assertIn("P3147/S2097", MD.read_text(encoding="utf-8"))
        self.assertIn("P3147/S2097", (ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md").read_text(encoding="utf-8"))
        self.assertIn("P3147/S2097", (ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md").read_text(encoding="utf-8"))
        self.assertIn("P3147", (REPO / "AGENTS.md").read_text(encoding="utf-8"))


if __name__ == "__main__":
    unittest.main()
