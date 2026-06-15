import json
import subprocess
import unittest
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]
ROOT = REPO / "fundamental_action_reconstruction"
SCRIPT = ROOT / "p2760_s1710_foundation_kernel_lagrangian_gap_matrix.py"
OUT = ROOT / "generated" / "p2760_s1710_foundation_kernel_lagrangian_gap_matrix.json"
MD = ROOT / "generated" / "p2760_s1710_foundation_kernel_lagrangian_gap_matrix.md"


class P2760FoundationKernelLagrangianGapMatrixTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        subprocess.run(["python", str(SCRIPT)], cwd=REPO, check=True)
        cls.payload = json.loads(OUT.read_text(encoding="utf-8"))

    def test_status_and_inputs(self):
        self.assertEqual(self.payload["status"], "P2760_FOUNDATION_KERNEL_LAGRANGIAN_GAP_MATRIX_NO_CLOSURE")
        self.assertEqual(self.payload["input_statuses"]["P1562"], "PASS_STRICT_KERNEL_COEFFICIENTS_DERIVED")
        self.assertEqual(self.payload["input_statuses"]["P1563"], "PASS_STRICT_KERNEL_TO_EOM_CHAIN_EXPORTED")
        self.assertEqual(self.payload["input_statuses"]["P1866"], "OPEN_OBSTRUCTION_WITH_TRACE")

    def test_content_scan_and_kernel_difference(self):
        self.assertTrue(self.payload["content_evidence_scan"]["all_patterns_have_hits"])
        comparison = self.payload["finite_kernel_comparison"]
        self.assertFalse(comparison["same_formula_on_sample_domain"])
        self.assertEqual(comparison["sample_domain"], list(range(13)))
        self.assertGreater(comparison["max_abs_delta"], 0.0)
        self.assertGreater(comparison["amplitude_ratio_at_d0_legacy_over_strict"], 1.0)

    def test_gap_matrix_has_professorial_obligations(self):
        gaps = self.payload["formula_gap_matrix"]
        self.assertEqual(gaps["gap_count"], 7)
        self.assertEqual(gaps["open_gap_count"], 7)
        self.assertEqual(gaps["closed_gap_count"], 0)
        self.assertTrue(gaps["stale_closure_flag_detected"])
        ids = {row["gap_id"] for row in gaps["rows"]}
        self.assertIn("G5_kernel_moments_to_physical_couplings", ids)
        self.assertIn("G7_closure_flag_consistency", ids)

    def test_acceptance_blocks_closure(self):
        matrix = self.payload["acceptance_matrix"]
        self.assertFalse(matrix["accepted_as_foundation_to_lagrangian_closure"])
        self.assertIn("bridge_theorem_exported", matrix["missing_criteria"])
        self.assertIn("selector_or_p2721_source_exported", matrix["missing_criteria"])
        self.assertIn("ltotal_reverse_closure_exported", matrix["missing_criteria"])
        flags = self.payload["decision"]["negative_export_flags"]
        self.assertTrue(all(value is False for value in flags.values()))

    def test_documentation_and_recommendation(self):
        self.assertIn("G5", self.payload["decision"]["next_honest_step"])
        self.assertIn("P2697-P2760", self.payload["decision"]["next_honest_step"])
        self.assertIn("P2760/S1710", MD.read_text(encoding="utf-8"))
        self.assertIn("P2760/S1710", (ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md").read_text(encoding="utf-8"))
        self.assertIn("P2760/S1710", (ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md").read_text(encoding="utf-8"))
        self.assertIn("P2760", (REPO / "AGENTS.md").read_text(encoding="utf-8"))


if __name__ == "__main__":
    unittest.main()
