from __future__ import annotations

import json
import subprocess
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent
SCRIPT = ROOT / "p2658_s1608_local_homogeneous_action_quantization_scale_clock_no_go.py"
OUT = ROOT / "generated" / "p2658_s1608_local_homogeneous_action_quantization_scale_clock_no_go.json"
MD = ROOT / "generated" / "p2658_s1608_local_homogeneous_action_quantization_scale_clock_no_go.md"


class P2658LocalHomogeneousActionQuantizationNoGoTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        subprocess.run([sys.executable, str(SCRIPT)], cwd=ROOT, check=True)
        cls.payload = json.loads(OUT.read_text(encoding="utf-8"))

    def test_content_first_grep_and_upstream_consistency(self) -> None:
        audit = self.payload["semantic_rg_antiduplication_audit"]
        self.assertIn("content-first", audit["mode"])
        self.assertIn("homogeneous_local_operator_content", audit["patterns"])
        self.assertIn("scale_clock_no_go_content", audit["patterns"])
        self.assertIn("action_quantization_content", audit["patterns"])
        self.assertTrue(all(data["count"] > 0 for data in audit["patterns"].values()))
        upstream = self.payload["upstream_consistency"]
        self.assertTrue(upstream["p2656_operator_covariance_verified"])
        self.assertTrue(upstream["p2656_absolute_operator_scale_not_selected"])
        self.assertTrue(upstream["p2657_integer_phase_condition_satisfied_all_scales"])
        self.assertTrue(upstream["p2657_unique_scale_not_selected"])
        self.assertTrue(upstream["p2657_next_step_allows_broader_no_go"])

    def test_homogeneous_no_go_verifies_covariance_and_phase_family(self) -> None:
        no_go = self.payload["homogeneous_quantization_no_go"]
        self.assertTrue(no_go["all_homogeneous_covariances_verified"])
        self.assertTrue(no_go["all_integer_phase_conditions_satisfied"])
        self.assertLess(no_go["max_covariance_error"], 1e-11)
        self.assertLess(no_go["max_tau_ratio_error_from_expected_degree"], 1e-11)
        self.assertFalse(no_go["unique_scale_selected_by_local_homogeneous_quantization"])
        self.assertEqual(no_go["audited_functional_count"], 10)

    def test_fixed_clock_generalization_is_external(self) -> None:
        fixed_clock = self.payload["fixed_clock_selector_generalization"]
        self.assertEqual(fixed_clock["false_selector_count"], 10)
        self.assertTrue(fixed_clock["all_audited_fixed_clock_selectors_are_external"])
        self.assertFalse(fixed_clock["fixed_clock_selector_admissible_now"])

    def test_source_matrix_blocks_current_class_but_keeps_anomaly_target_open(self) -> None:
        matrix = {row["candidate"]: row for row in self.payload["source_candidate_matrix"]}
        self.assertFalse(matrix["local_homogeneous_integer_phase_family"]["passes_as_uv_unit_source_now"])
        self.assertFalse(matrix["local_homogeneous_integer_phase_family"]["selects_unique_scale"])
        self.assertTrue(matrix["fixed_clock_for_homogeneous_functional"]["uses_external_clock_or_scale_anchor"])
        self.assertFalse(matrix["fixed_clock_for_homogeneous_functional"]["passes_as_uv_unit_source_now"])
        self.assertFalse(matrix["nonhomogeneous_or_boundary_anomaly_source"]["covered_by_audit"])
        self.assertFalse(matrix["intrinsic_nadsoliton_clock_source_theorem"]["passes_as_uv_unit_source_now"])

    def test_no_closure_and_docs_updated(self) -> None:
        decision = self.payload["closure_decision"]
        self.assertEqual(decision["passing_uv_unit_source_candidates"], [])
        self.assertTrue(decision["all_homogeneous_covariances_verified"])
        self.assertTrue(decision["all_integer_phase_conditions_satisfied"])
        self.assertFalse(decision["uv_unit_selected_now"])
        self.assertFalse(decision["beta_source_exported_now"])
        self.assertFalse(decision["role_bearing_ltotal_now"])
        self.assertFalse(decision["toe_closure_now"])
        self.assertTrue(all(value is False for value in self.payload["negative_export_flags"].values()))
        self.assertIn("P2658/S1608", MD.read_text(encoding="utf-8"))
        self.assertIn("P2658/S1608", (ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md").read_text(encoding="utf-8"))
        self.assertIn("P2658/S1608", (ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md").read_text(encoding="utf-8"))


if __name__ == "__main__":
    unittest.main()
