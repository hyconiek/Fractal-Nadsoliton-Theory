from __future__ import annotations

import json
import subprocess
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent
SCRIPT = ROOT / "p2657_s1607_nadsoliton_action_quantization_scale_anchor_obstruction_audit.py"
OUT = ROOT / "generated" / "p2657_s1607_nadsoliton_action_quantization_scale_anchor_obstruction_audit.json"
MD = ROOT / "generated" / "p2657_s1607_nadsoliton_action_quantization_scale_anchor_obstruction_audit.md"


class P2657NadsolitonActionQuantizationScaleAnchorAuditTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        subprocess.run([sys.executable, str(SCRIPT)], cwd=ROOT, check=True)
        cls.payload = json.loads(OUT.read_text(encoding="utf-8"))

    def test_content_first_grep_and_upstream_consistency(self) -> None:
        audit = self.payload["semantic_rg_antiduplication_audit"]
        self.assertIn("content-first", audit["mode"])
        self.assertIn("action_quantization_content", audit["patterns"])
        self.assertIn("clock_scale_anchor_content", audit["patterns"])
        self.assertIn("uv_beta_source_content", audit["patterns"])
        self.assertTrue(all(data["count"] > 0 for data in audit["patterns"].values()))
        upstream = self.payload["upstream_consistency"]
        self.assertTrue(upstream["p2655_uv_unit_not_selected"])
        self.assertTrue(upstream["p2655_scale_quotient_rank_one"])
        self.assertTrue(upstream["p2656_operator_covariance_verified"])
        self.assertTrue(upstream["p2656_absolute_operator_scale_not_selected"])
        self.assertTrue(upstream["p2656_next_step_requests_action_quantization"])

    def test_integer_phase_condition_leaves_scale_clock_family(self) -> None:
        family = self.payload["action_quantization_family_audit"]
        self.assertTrue(family["integer_phase_condition_satisfied_all_scales"])
        self.assertLess(family["max_trace_covariance_error"], 1e-12)
        self.assertLess(family["max_tau_ratio_error_from_a_squared"], 1e-12)
        self.assertFalse(family["unique_scale_selected_by_integer_phase_alone"])
        self.assertGreater(family["quantized_family_count"], 1)

    def test_fixed_clock_false_selector_is_blocked(self) -> None:
        fixed_clock = self.payload["fixed_clock_false_selector_audit"]
        self.assertTrue(fixed_clock["scale_one_only_passes_with_imported_clock"])
        self.assertTrue(fixed_clock["fixed_clock_is_external_anchor"])
        self.assertFalse(fixed_clock["fixed_clock_selector_admissible_now"])

    def test_source_matrix_blocks_current_quantization_candidates(self) -> None:
        matrix = {row["candidate"]: row for row in self.payload["source_candidate_matrix"]}
        self.assertFalse(matrix["integer_phase_condition_tau_trace_l_equals_2pi_n"]["selects_unique_scale"])
        self.assertFalse(matrix["integer_phase_condition_tau_trace_l_equals_2pi_n"]["passes_as_uv_unit_source_now"])
        self.assertTrue(matrix["fixed_tau_from_scale_one_selector"]["uses_external_clock_or_scale_anchor"])
        self.assertFalse(matrix["fixed_tau_from_scale_one_selector"]["passes_as_uv_unit_source_now"])
        self.assertTrue(matrix["declared_trace_or_gap_action_quantum"]["uses_external_clock_or_scale_anchor"])
        self.assertFalse(matrix["dimensionless_integer_phase_ratio_identity"]["selects_unique_scale"])
        self.assertFalse(matrix["full_intrinsic_nadsoliton_quantization_theorem"]["passes_as_uv_unit_source_now"])

    def test_no_closure_and_docs_updated(self) -> None:
        decision = self.payload["closure_decision"]
        self.assertEqual(decision["passing_uv_unit_source_candidates"], [])
        self.assertTrue(decision["integer_phase_condition_satisfied_all_scales"])
        self.assertFalse(decision["unique_scale_selected_by_integer_phase_alone"])
        self.assertFalse(decision["fixed_clock_selector_admissible_now"])
        self.assertFalse(decision["uv_unit_selected_now"])
        self.assertFalse(decision["beta_source_exported_now"])
        self.assertFalse(decision["role_bearing_ltotal_now"])
        self.assertFalse(decision["toe_closure_now"])
        self.assertTrue(all(value is False for value in self.payload["negative_export_flags"].values()))
        self.assertIn("P2657/S1607", MD.read_text(encoding="utf-8"))
        self.assertIn("P2657/S1607", (ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md").read_text(encoding="utf-8"))
        self.assertIn("P2657/S1607", (ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md").read_text(encoding="utf-8"))


if __name__ == "__main__":
    unittest.main()
