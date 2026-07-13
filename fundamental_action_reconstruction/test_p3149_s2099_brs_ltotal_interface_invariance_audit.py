import json
import subprocess
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
SCRIPT = ROOT / "p3149_s2099_brs_ltotal_interface_invariance_audit.py"
OUT = ROOT / "generated" / "p3149_s2099_brs_ltotal_interface_invariance_audit.json"
MD = ROOT / "generated" / "p3149_s2099_brs_ltotal_interface_invariance_audit.md"


class P3149BrstLtotalInterfaceInvarianceAuditTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        subprocess.run([sys.executable, str(SCRIPT)], cwd=ROOT, check=True, stdout=subprocess.PIPE, text=True)
        cls.payload = json.loads(OUT.read_text(encoding="utf-8"))

    def test_status_and_counts(self):
        self.assertEqual(self.payload["status"], "P3149_BRST_LTOTAL_INTERFACE_LOCAL_ALGEBRAIC_PASS_CONDITIONAL_NO_GLOBAL_CLOSURE")
        counts = self.payload["finite_theorem"]["finite_counts"]
        self.assertEqual(counts["kinetic_rows"], 6)
        self.assertEqual(counts["kinetic_rows_local_gauge_invariant"], 6)
        self.assertEqual(counts["yukawa_rows"], 3)
        self.assertEqual(counts["yukawa_rows_local_brst_zero"], 3)
        self.assertEqual(counts["local_lagrangian_interface_rows"], 9)
        self.assertEqual(counts["local_lagrangian_interface_rows_passing"], 9)
        self.assertEqual(counts["unit_bearing_measure_rows"], 0)
        self.assertEqual(counts["global_bv_brst_rows"], 0)
        self.assertEqual(counts["strict_nadsoliton_source_rows"], 0)

    def test_invariance_rows_pass(self):
        self.assertTrue(all(row["local_gauge_invariant_bilinear"] for row in self.payload["kinetic_invariance_rows"]))
        self.assertTrue(all(row["local_brst_variation_zero"] for row in self.payload["yukawa_invariance_rows"]))
        terms = {row["term"] for row in self.payload["yukawa_invariance_rows"]}
        self.assertEqual(terms, {"Q_L H u_c", "Q_L Hdagger d_c", "L_L Hdagger e_c"})

    def test_decision_preserves_global_boundaries(self):
        decision = self.payload["decision"]
        self.assertTrue(decision["accepted_as_local_interface_certificate"])
        self.assertIn("local and algebraic", decision["why_not_strict"])
        self.assertIn("P3150", decision["next_honest_step"])
        self.assertTrue(all(value is False for value in decision["negative_export_flags"].values()))

    def test_docs_updated(self):
        self.assertIn("P3149/S2099", MD.read_text(encoding="utf-8"))
        self.assertIn("P3149/S2099", (ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md").read_text(encoding="utf-8"))
        self.assertIn("P3149/S2099", (ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md").read_text(encoding="utf-8"))
        self.assertIn("P3149", (REPO / "AGENTS.md").read_text(encoding="utf-8"))


if __name__ == "__main__":
    unittest.main()
