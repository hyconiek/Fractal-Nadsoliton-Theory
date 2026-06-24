import json, subprocess, sys, unittest
from pathlib import Path
ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
SCRIPT = ROOT / "p3072_s2022_discrete_continuity_noether_current_matrix.py"
OUT = ROOT / "generated" / "p3072_s2022_discrete_continuity_noether_current_matrix.json"
MD = ROOT / "generated" / "p3072_s2022_discrete_continuity_noether_current_matrix.md"

class P3072DiscreteContinuityNoetherCurrentMatrixTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        subprocess.run([sys.executable, str(SCRIPT)], cwd=ROOT, check=True, stdout=subprocess.PIPE, text=True)
        cls.payload = json.loads(OUT.read_text(encoding="utf-8"))

    def test_status_and_input_hash(self):
        self.assertEqual(self.payload["status"], "P3072_DISCRETE_CONTINUITY_NOETHER_CURRENT_INTERFACE_BOUNDED_NO_GO")
        self.assertIsNotNone(self.payload["input_hashes"]["P3071"])

    def test_content_first_and_matrix_counts(self):
        cert = self.payload["finite_certificate"]
        self.assertEqual(cert["content_grep_lanes"], 4)
        self.assertGreater(cert["content_grep_hits"], 0)
        self.assertEqual(cert["p3071_accepted_profiles"], 3)
        self.assertEqual(cert["profiles_tested"], 3)
        self.assertEqual(cert["sigma_branches"], 2)
        self.assertEqual(cert["d12_transforms"], 24)
        self.assertEqual(cert["current_templates"], 5)
        self.assertEqual(cert["continuity_matrix_rows"], 720)
        self.assertEqual(cert["divergence_zero_rows"], 480)
        self.assertEqual(cert["orientation_premise_rows"], 288)
        self.assertEqual(cert["accepted_premise_free_static_continuity_rows"], 144)
        self.assertEqual(cert["accepted_nontrivial_noether_current_rows"], 0)
        self.assertEqual(cert["nonzero_divergence_rows"], 240)
        self.assertEqual(cert["satisfied_proof_obligations"], 4)

    def test_current_templates_and_obstruction(self):
        rows = self.payload["constructed_theoretical_objects"]["continuity_current_matrix"]
        accepted = [r for r in rows if r["accepted_premise_free_static_continuity_row"]]
        nontrivial = [r for r in rows if r["accepted_nontrivial_noether_current_row"]]
        oriented = [r for r in rows if r["orientation_premise_required"] and r["divergence_zero"]]
        residual = [r for r in rows if not r["divergence_zero"]]
        self.assertEqual(len(accepted), 144)
        self.assertTrue(all(r["current_template"] == "zero_current" for r in accepted))
        self.assertEqual(len(nontrivial), 0)
        self.assertEqual(len(oriented), 288)
        self.assertEqual(len(residual), 240)

    def test_flags_recommendation_and_docs(self):
        for value in self.payload["decision"]["negative_export_flags"].values():
            self.assertFalse(value)
        self.assertTrue(self.payload["decision"]["positive_scoped_flags"]["z12_incidence_divergence_operator_exported"])
        self.assertIn("bounded renormalization/scale-flow obstruction table", self.payload["decision"]["next_honest_step"])
        self.assertIn("P3072/S2022", MD.read_text(encoding="utf-8"))
        self.assertIn("P3072/S2022", (ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md").read_text(encoding="utf-8"))
        self.assertIn("P3072/S2022", (ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md").read_text(encoding="utf-8"))
        self.assertIn("P3072", (REPO / "AGENTS.md").read_text(encoding="utf-8"))

if __name__ == "__main__":
    unittest.main()
