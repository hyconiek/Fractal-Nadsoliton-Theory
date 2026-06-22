import json, subprocess, sys, unittest
from pathlib import Path
ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
SCRIPT = ROOT / "p3012_s1962_strict_kernel_phase_gradient_selector_source_obstruction.py"
OUT = ROOT / "generated" / "p3012_s1962_strict_kernel_phase_gradient_selector_source_obstruction.json"
MD = ROOT / "generated" / "p3012_s1962_strict_kernel_phase_gradient_selector_source_obstruction.md"

class P3012StrictKernelPhaseGradientSelectorSourceTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        subprocess.run([sys.executable, str(SCRIPT)], cwd=ROOT, check=True, stdout=subprocess.PIPE, text=True)
        cls.payload = json.loads(OUT.read_text(encoding="utf-8"))

    def test_status_and_input(self):
        self.assertEqual(self.payload["status"], "P3012_STRICT_KERNEL_PHASE_GRADIENT_SELECTOR_SOURCE_OBSTRUCTION_BOUNDED_NO_GO")
        self.assertIsNotNone(self.payload["input_hashes"]["P3011"])

    def test_selector_certificate(self):
        cert = self.payload["selector_certificate"]
        self.assertEqual(cert["label_based_pair_pick"], 1)
        self.assertIn(cert["label_based_unit_pick"], [1, 5, 7, 11])
        self.assertNotEqual(cert["pair_score_gap_K1_minus_K11"], 0)
        self.assertEqual(cert["equivariance_row_count"], 16)
        self.assertGreater(cert["equivariance_failure_count"], 0)
        self.assertGreater(cert["pair_action_leaves_pair_count"], 0)
        self.assertEqual(cert["aut_invariant_directed_unit_count"], 0)
        self.assertEqual(cert["acceptance_matrix_rows"], 128)
        self.assertEqual(cert["accepted_rows"], 1)

    def test_obligations_and_nonpromotion(self):
        obj = self.payload["constructed_theoretical_objects"]
        self.assertEqual(obj["object"], "StrictKernelPhaseGradientSelectorSource_EquivarianceObstructionMatrix")
        obligations = {r["obligation"]: r["satisfied"] for r in obj["proof_obligation_rows"]}
        self.assertTrue(obligations["finite_label_arrow_exists"])
        self.assertTrue(obligations["unit_orbit_score_table_exists"])
        self.assertFalse(obligations["aut_equivariance"])
        self.assertFalse(obligations["pair_closed_under_full_aut_action"])
        self.assertFalse(obligations["nonpremise_chart_metric_source"])
        self.assertFalse(obligations["accepted_current_selector_source"])
        for value in self.payload["decision"]["negative_export_flags"].values():
            self.assertFalse(value)

    def test_docs(self):
        self.assertIn("P3012/S1962", MD.read_text(encoding="utf-8"))
        self.assertIn("P3012/S1962", (ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md").read_text(encoding="utf-8"))
        self.assertIn("P3012/S1962", (ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md").read_text(encoding="utf-8"))
        self.assertIn("P3012", (REPO / "AGENTS.md").read_text(encoding="utf-8"))

if __name__ == "__main__":
    unittest.main()
