import json, subprocess, sys, unittest
from pathlib import Path
ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
SCRIPT = ROOT / "p3065_s2015_fundamental_selector_axiom_boundary_certificate.py"
OUT = ROOT / "generated" / "p3065_s2015_fundamental_selector_axiom_boundary_certificate.json"
MD = ROOT / "generated" / "p3065_s2015_fundamental_selector_axiom_boundary_certificate.md"

class P3065FundamentalSelectorAxiomBoundaryTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        subprocess.run([sys.executable, str(SCRIPT)], cwd=ROOT, check=True, stdout=subprocess.PIPE, text=True)
        cls.payload = json.loads(OUT.read_text(encoding="utf-8"))

    def test_status_and_input(self):
        self.assertEqual(self.payload["status"], "P3065_FUNDAMENTAL_SELECTOR_AXIOM_BOUNDARY_CERTIFICATE_AXIOM_AUGMENTED_ONLY")
        self.assertIsNotNone(self.payload["input_hashes"]["P3064"])

    def test_content_first_and_boundary_counts(self):
        cert = self.payload["finite_certificate"]
        self.assertEqual(cert["content_grep_lanes"], 3)
        self.assertGreater(cert["content_grep_hits"], 0)
        self.assertEqual(cert["sigma_choices"], 2)
        self.assertEqual(cert["boundary_rows"], 3)
        self.assertEqual(cert["axiom_augmented_admitted_rows"], 2)
        self.assertEqual(cert["strict_selector_export_rows"], 0)
        self.assertEqual(cert["current_strict_exported_atoms"], 0)
        self.assertEqual(cert["conditioned_nadsoliton_continuation_rows"], 2)
        self.assertFalse(cert["p3064_current_row_accepted"])
        self.assertEqual(cert["satisfied_proof_obligations"], 4)

    def test_constructed_axiom_object_and_rows(self):
        obj = self.payload["constructed_theoretical_objects"]["fundamental_selector_axiom"]
        self.assertEqual(obj["object"], "FundamentalSelectorAxiomBoundary")
        self.assertIn("boundary parameter", obj["interpretation"])
        rows = self.payload["constructed_theoretical_objects"]["axiom_boundary_rows"]
        admitted = [row for row in rows if row["axiom_augmented_theory_admitted"]]
        self.assertEqual({row["sigma_axiom"] for row in admitted}, {"sigma_plus", "sigma_minus"})
        self.assertTrue(all(row["may_continue_nadsoliton_analysis_conditioned_on_sigma"] for row in admitted))
        self.assertTrue(all(not row["strict_selector_export"] for row in rows))

    def test_export_flags_and_docs(self):
        for value in self.payload["decision"]["negative_export_flags"].values():
            self.assertFalse(value)
        self.assertTrue(self.payload["decision"]["positive_scoped_flags"]["axiom_augmented_T_sigma_family_constructed"])
        self.assertIn("P3065/S2015", MD.read_text(encoding="utf-8"))
        self.assertIn("P3065/S2015", (ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md").read_text(encoding="utf-8"))
        self.assertIn("P3065/S2015", (ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md").read_text(encoding="utf-8"))
        self.assertIn("P3065", (REPO / "AGENTS.md").read_text(encoding="utf-8"))

if __name__ == "__main__":
    unittest.main()
