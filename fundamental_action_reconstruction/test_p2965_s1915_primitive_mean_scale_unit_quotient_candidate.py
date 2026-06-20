import json, subprocess, sys, unittest
from pathlib import Path
ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
SCRIPT = ROOT / "p2965_s1915_primitive_mean_scale_unit_quotient_candidate.py"
OUT = ROOT / "generated" / "p2965_s1915_primitive_mean_scale_unit_quotient_candidate.json"
MD = ROOT / "generated" / "p2965_s1915_primitive_mean_scale_unit_quotient_candidate.md"

class P2965PrimitiveMeanScaleUnitQuotientCandidateTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        subprocess.run([sys.executable, str(SCRIPT)], cwd=ROOT, check=True, stdout=subprocess.PIPE, text=True)
        cls.payload = json.loads(OUT.read_text(encoding="utf-8"))

    def test_status_and_input(self):
        self.assertEqual(self.payload["status"], "P2965_PRIMITIVE_MEAN_SCALE_UNIT_QUOTIENT_CANDIDATE_NO_STRICT_EXPORT")
        self.assertIsNotNone(self.payload["input_hashes"]["P2964"])
        self.assertIn("primitive mean", self.payload["unbounded_lemma"])

    def test_ray_samples_and_candidates(self):
        cert = self.payload["quotient_certificate"]
        self.assertEqual(cert["sample_count"], 7)
        self.assertTrue(cert["all_samples_invariant"])
        self.assertEqual(cert["developmental_candidates"], ["primitive_ray_mean_quotient"])
        self.assertEqual(cert["strict_coupling_sources"], [])
        for row in self.payload["constructed_theoretical_objects"]["ray_sample_rows"]:
            self.assertEqual(row["primitive_representative"], [1, 2, 2, 2, 2])
            self.assertEqual(row["primitive_sum"], 9)
            self.assertEqual(row["primitive_mean"], "9/5")

    def test_obligations_and_matrix(self):
        obligations = {r["obligation"]: r["satisfied"] for r in self.payload["constructed_theoretical_objects"]["proof_obligation_rows"]}
        self.assertTrue(obligations["positive_ray_quotient_constructed"])
        self.assertTrue(obligations["finite_scale_samples_invariant"])
        self.assertTrue(obligations["coefficient_9_over_5_computed_without_sum_cut"])
        self.assertFalse(obligations["physical_unit_law_exported"])
        self.assertFalse(obligations["strict_coefficient_source_law_exported"])
        self.assertFalse(obligations["nonproxy_coupling_installed"])
        matrix = self.payload["constructed_theoretical_objects"]["finite_acceptance_matrix"]
        self.assertEqual(len(matrix), 128)
        self.assertEqual(sum(1 for row in matrix if row["accepts_developmental_quotient_candidate"]), 8)
        self.assertEqual(sum(1 for row in matrix if row["accepts_strict_coupling_source"]), 1)

    def test_docs_and_nonpromotion(self):
        for value in self.payload["decision"]["negative_export_flags"].values():
            self.assertFalse(value)
        self.assertIn("P2965/S1915", MD.read_text(encoding="utf-8"))
        self.assertIn("P2965/S1915", (ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md").read_text(encoding="utf-8"))
        self.assertIn("P2965/S1915", (ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md").read_text(encoding="utf-8"))
        self.assertIn("P2965", (REPO / "AGENTS.md").read_text(encoding="utf-8"))

if __name__ == "__main__":
    unittest.main()
