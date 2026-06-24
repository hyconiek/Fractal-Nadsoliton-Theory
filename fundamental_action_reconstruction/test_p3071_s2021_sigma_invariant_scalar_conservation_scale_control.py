import json, subprocess, sys, unittest
from pathlib import Path
ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
SCRIPT = ROOT / "p3071_s2021_sigma_invariant_scalar_conservation_scale_control.py"
OUT = ROOT / "generated" / "p3071_s2021_sigma_invariant_scalar_conservation_scale_control.json"
MD = ROOT / "generated" / "p3071_s2021_sigma_invariant_scalar_conservation_scale_control.md"

class P3071SigmaInvariantScalarConservationScaleControlTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        subprocess.run([sys.executable, str(SCRIPT)], cwd=ROOT, check=True, stdout=subprocess.PIPE, text=True)
        cls.payload = json.loads(OUT.read_text(encoding="utf-8"))

    def test_status_and_input_hash(self):
        self.assertEqual(self.payload["status"], "P3071_SIGMA_INVARIANT_SCALAR_CONSERVATION_SCALE_CONTROL_SCOPED_EXPORT")
        self.assertIsNotNone(self.payload["input_hashes"]["P3070"])

    def test_content_first_and_matrix_counts(self):
        cert = self.payload["finite_certificate"]
        self.assertEqual(cert["content_grep_lanes"], 4)
        self.assertGreater(cert["content_grep_hits"], 0)
        self.assertEqual(cert["sigma_branches"], 2)
        self.assertEqual(cert["z12_rows"], 12)
        self.assertEqual(cert["d12_transforms"], 24)
        self.assertEqual(cert["candidate_profiles"], 5)
        self.assertEqual(cert["conservation_matrix_rows"], 240)
        self.assertEqual(cert["accepted_profile_count"], 3)
        self.assertEqual(cert["accepted_matrix_rows"], 144)
        self.assertEqual(cert["rejected_matrix_rows"], 96)
        self.assertEqual(cert["p3070_accepted_unit_provider_rows"], 0)
        self.assertEqual(cert["satisfied_proof_obligations"], 4)

    def test_accepted_and_rejected_profiles(self):
        cert = self.payload["finite_certificate"]
        self.assertEqual(cert["accepted_profile_ids"], [
            "constant_cardinality_density",
            "even_distance_quadratic_density",
            "even_distance_shell_indicator_density",
        ])
        aggregates = self.payload["constructed_theoretical_objects"]["profile_aggregate_certificate"]
        accepted = [row for row in aggregates if row["all_rows_accepted"]]
        rejected = [row for row in aggregates if not row["all_rows_accepted"]]
        self.assertEqual(len(accepted), 3)
        self.assertEqual(len(rejected), 2)
        self.assertTrue(all(not row["exports_observed_physics"] for row in aggregates))

    def test_flags_recommendation_and_docs(self):
        for value in self.payload["decision"]["negative_export_flags"].values():
            self.assertFalse(value)
        self.assertTrue(self.payload["decision"]["positive_scoped_flags"]["bounded_dimensionless_scale_control_exported"])
        self.assertIn("discrete continuity/Noether-current matrix", self.payload["decision"]["next_honest_step"])
        self.assertIn("P3071/S2021", MD.read_text(encoding="utf-8"))
        self.assertIn("P3071/S2021", (ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md").read_text(encoding="utf-8"))
        self.assertIn("P3071/S2021", (ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md").read_text(encoding="utf-8"))
        self.assertIn("P3071", (REPO / "AGENTS.md").read_text(encoding="utf-8"))

if __name__ == "__main__":
    unittest.main()
