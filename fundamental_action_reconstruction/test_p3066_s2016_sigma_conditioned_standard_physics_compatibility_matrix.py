import json, subprocess, sys, unittest
from pathlib import Path
ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
SCRIPT = ROOT / "p3066_s2016_sigma_conditioned_standard_physics_compatibility_matrix.py"
OUT = ROOT / "generated" / "p3066_s2016_sigma_conditioned_standard_physics_compatibility_matrix.json"
MD = ROOT / "generated" / "p3066_s2016_sigma_conditioned_standard_physics_compatibility_matrix.md"

class P3066SigmaConditionedStandardPhysicsCompatibilityTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        subprocess.run([sys.executable, str(SCRIPT)], cwd=ROOT, check=True, stdout=subprocess.PIPE, text=True)
        cls.payload = json.loads(OUT.read_text(encoding="utf-8"))

    def test_status_and_input_hash(self):
        self.assertEqual(self.payload["status"], "P3066_SIGMA_CONDITIONED_STANDARD_PHYSICS_COMPATIBILITY_MATRIX_NO_CLOSURE")
        self.assertIsNotNone(self.payload["input_hashes"]["P3065"])

    def test_content_first_and_finite_counts(self):
        cert = self.payload["finite_certificate"]
        self.assertEqual(cert["content_grep_lanes"], 3)
        self.assertGreater(cert["content_grep_hits"], 0)
        self.assertEqual(cert["sigma_branches"], 2)
        self.assertEqual(cert["observable_rows"], 6)
        self.assertEqual(cert["branch_matrix_rows"], 12)
        self.assertEqual(cert["sigma_invariant_observables"], 2)
        self.assertEqual(cert["sigma_odd_observables"], 2)
        self.assertEqual(cert["unknown_parity_observables"], 2)
        self.assertEqual(cert["physics_obligations"], 6)
        self.assertEqual(cert["physics_obligation_rows"], 36)
        self.assertEqual(cert["accepted_physics_obligation_rows"], 0)
        self.assertEqual(cert["satisfied_proof_obligations"], 4)

    def test_constructed_matrix(self):
        matrix = self.payload["constructed_theoretical_objects"]["sigma_conditioned_standard_physics_compatibility_matrix"]
        self.assertEqual(matrix["object"], "SigmaConditionedStandardPhysicsCompatibilityMatrix")
        self.assertEqual(len(matrix["sigma_branches"]), 2)
        self.assertEqual(len(matrix["physics_obligations"]), 6)
        self.assertEqual(len(matrix["observable_rows"]), 6)
        self.assertEqual(len(matrix["branch_matrix"]), 12)
        self.assertEqual(len(matrix["physics_obligation_matrix"]), 36)
        odd_rows = [row for row in matrix["observable_rows"] if row["sigma_parity"] == "odd"]
        self.assertTrue(all("requires_explicit_pseudoscalar_or_orientation_coupling_for_standard_physics_use" in row["remaining_blockers"] for row in odd_rows))

    def test_flags_recommendation_and_docs(self):
        for value in self.payload["decision"]["negative_export_flags"].values():
            self.assertFalse(value)
        self.assertTrue(self.payload["decision"]["positive_scoped_flags"]["selector_search_no_longer_blocks_conditioned_analysis"])
        self.assertIn("light_emergence_interface", self.payload["decision"]["next_honest_step"])
        self.assertIn("P3066/S2016", MD.read_text(encoding="utf-8"))
        self.assertIn("P3066/S2016", (ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md").read_text(encoding="utf-8"))
        self.assertIn("P3066/S2016", (ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md").read_text(encoding="utf-8"))
        self.assertIn("P3066", (REPO / "AGENTS.md").read_text(encoding="utf-8"))

if __name__ == "__main__":
    unittest.main()
