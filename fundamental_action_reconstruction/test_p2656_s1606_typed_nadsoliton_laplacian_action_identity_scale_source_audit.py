from __future__ import annotations

import json
import subprocess
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent
SCRIPT = ROOT / "p2656_s1606_typed_nadsoliton_laplacian_action_identity_scale_source_audit.py"
OUT = ROOT / "generated" / "p2656_s1606_typed_nadsoliton_laplacian_action_identity_scale_source_audit.json"
MD = ROOT / "generated" / "p2656_s1606_typed_nadsoliton_laplacian_action_identity_scale_source_audit.md"


class P2656TypedLaplacianActionIdentityScaleSourceAuditTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        subprocess.run([sys.executable, str(SCRIPT)], cwd=ROOT, check=True)
        cls.payload = json.loads(OUT.read_text(encoding="utf-8"))

    def test_content_first_grep_and_upstream_consistency(self) -> None:
        audit = self.payload["semantic_rg_antiduplication_audit"]
        self.assertIn("content-first", audit["mode"])
        self.assertIn("typed_dynamics_operator_content", audit["patterns"])
        self.assertIn("scale_source_content", audit["patterns"])
        self.assertIn("beta_source_content", audit["patterns"])
        self.assertTrue(all(data["count"] > 0 for data in audit["patterns"].values()))
        upstream = self.payload["upstream_consistency"]
        self.assertTrue(upstream["p2653_typed_metric_uv_source_not_exported"])
        self.assertTrue(upstream["p2654_scale_invariant_selector_blocked"])
        self.assertTrue(upstream["p2655_metric_skeleton_constructed"])
        self.assertTrue(upstream["p2655_uv_unit_not_selected"])
        self.assertTrue(upstream["p2655_next_step_requests_operator_identity"])

    def test_laplacian_operator_covariant_but_dimensionless_rank_one(self) -> None:
        operator = self.payload["laplacian_scale_covariance_audit"]
        self.assertTrue(operator["operator_covariance_verified"])
        self.assertLess(operator["max_covariance_error"], 1e-12)
        self.assertEqual(operator["dimensionless_operator_rank_on_scale_orbit"], 1)
        self.assertLess(operator["max_dimensionless_invariant_error"], 1e-12)
        self.assertFalse(operator["absolute_operator_scale_selected"])

    def test_operator_beta_covariance_does_not_source_beta_one(self) -> None:
        beta_audit = self.payload["beta_covariance_from_operator_scale_audit"]
        self.assertTrue(beta_audit["all_nonunit_scales_keep_beta_one_false"])
        self.assertTrue(beta_audit["trace_anchor_is_external_until_dynamics_quantizes_it"])
        self.assertFalse(beta_audit["beta_source_exported_by_operator_now"])

    def test_source_matrix_blocks_absolute_anchor_promotions(self) -> None:
        matrix = {row["candidate"]: row for row in self.payload["source_candidate_matrix"]}
        self.assertFalse(matrix["dimensionless_laplacian_trace_moment_identity"]["breaks_scale_orbit"])
        self.assertFalse(matrix["dimensionless_laplacian_trace_moment_identity"]["passes_as_uv_unit_source_now"])
        self.assertTrue(matrix["absolute_trace_l_equals_declared_quantum"]["uses_external_absolute_anchor"])
        self.assertFalse(matrix["absolute_trace_l_equals_declared_quantum"]["passes_as_uv_unit_source_now"])
        self.assertTrue(matrix["spectral_gap_or_heat_time_normalization"]["uses_external_absolute_anchor"])
        self.assertFalse(matrix["operator_beta_covariance_identity"]["passes_as_uv_unit_source_now"])
        self.assertFalse(matrix["full_nadsoliton_action_quantization_theorem"]["passes_as_uv_unit_source_now"])

    def test_no_closure_and_docs_updated(self) -> None:
        decision = self.payload["closure_decision"]
        self.assertEqual(decision["passing_uv_unit_source_candidates"], [])
        self.assertTrue(decision["operator_covariance_verified"])
        self.assertTrue(decision["dimensionless_operator_rank_one"])
        self.assertFalse(decision["absolute_operator_scale_selected"])
        self.assertFalse(decision["uv_unit_selected_now"])
        self.assertFalse(decision["beta_source_exported_now"])
        self.assertFalse(decision["role_bearing_ltotal_now"])
        self.assertFalse(decision["toe_closure_now"])
        self.assertTrue(all(value is False for value in self.payload["negative_export_flags"].values()))
        self.assertIn("P2656/S1606", MD.read_text(encoding="utf-8"))
        self.assertIn("P2656/S1606", (ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md").read_text(encoding="utf-8"))
        self.assertIn("P2656/S1606", (ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md").read_text(encoding="utf-8"))


if __name__ == "__main__":
    unittest.main()
