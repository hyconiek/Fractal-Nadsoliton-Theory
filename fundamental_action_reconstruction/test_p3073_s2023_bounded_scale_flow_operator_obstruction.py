import json, subprocess, sys, unittest
from pathlib import Path
ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
SCRIPT = ROOT / "p3073_s2023_bounded_scale_flow_operator_obstruction.py"
OUT = ROOT / "generated" / "p3073_s2023_bounded_scale_flow_operator_obstruction.json"
MD = ROOT / "generated" / "p3073_s2023_bounded_scale_flow_operator_obstruction.md"

class P3073BoundedScaleFlowOperatorObstructionTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        subprocess.run([sys.executable, str(SCRIPT)], cwd=ROOT, check=True, stdout=subprocess.PIPE, text=True)
        cls.payload = json.loads(OUT.read_text(encoding="utf-8"))

    def test_status_and_input_hash(self):
        self.assertEqual(self.payload["status"], "P3073_BOUNDED_SCALE_FLOW_OPERATOR_PARTIAL_EXPORT_FULL_SUMMARY_DYNAMICS_NO_GO")
        self.assertIsNotNone(self.payload["input_hashes"]["P3072"])

    def test_matrix_counts(self):
        cert = self.payload["finite_certificate"]
        self.assertEqual(cert["content_grep_lanes"], 4)
        self.assertGreater(cert["content_grep_hits"], 0)
        self.assertEqual(cert["p3072_accepted_nontrivial_noether_current_rows"], 0)
        self.assertEqual(cert["profiles_tested"], 3)
        self.assertEqual(cert["sigma_branches"], 2)
        self.assertEqual(cert["d12_transforms"], 24)
        self.assertEqual(cert["flow_operators"], 5)
        self.assertEqual(cert["scale_flow_matrix_rows"], 720)
        self.assertEqual(cert["accepted_intrinsic_bounded_scale_flow_rows"], 192)
        self.assertEqual(cert["accepted_full_summary_dynamics_rows"], 0)
        self.assertEqual(cert["total_preserving_rows"], 432)
        self.assertEqual(cert["nonzero_flow_rows"], 480)
        self.assertEqual(cert["satisfied_proof_obligations"], 4)

    def test_operator_aggregates(self):
        aggregates = self.payload["constructed_theoretical_objects"]["operator_aggregate_certificate"]
        by_id = {row["flow_operator"]: row for row in aggregates}
        self.assertEqual(by_id["cycle_laplacian_flow"]["accepted_intrinsic_bounded_scale_flow_rows"], 96)
        self.assertEqual(by_id["mean_centering_flow"]["accepted_intrinsic_bounded_scale_flow_rows"], 96)
        self.assertEqual(by_id["zero_flow"]["accepted_intrinsic_bounded_scale_flow_rows"], 0)
        self.assertEqual(sum(row["accepted_full_summary_dynamics_rows"] for row in aggregates), 0)

    def test_flags_recommendation_and_docs(self):
        for value in self.payload["decision"]["negative_export_flags"].values():
            self.assertFalse(value)
        self.assertTrue(self.payload["decision"]["positive_scoped_flags"]["intrinsic_bounded_total_preserving_scale_flow_exported"])
        self.assertIn("Lyapunov/entropy monotonicity certificate", self.payload["decision"]["next_honest_step"])
        self.assertIn("P3073/S2023", MD.read_text(encoding="utf-8"))
        self.assertIn("P3073/S2023", (ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md").read_text(encoding="utf-8"))
        self.assertIn("P3073/S2023", (ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md").read_text(encoding="utf-8"))
        self.assertIn("P3073", (REPO / "AGENTS.md").read_text(encoding="utf-8"))

if __name__ == "__main__":
    unittest.main()
