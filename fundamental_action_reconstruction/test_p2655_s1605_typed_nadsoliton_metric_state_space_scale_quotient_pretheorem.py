from __future__ import annotations

import json
import subprocess
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent
SCRIPT = ROOT / "p2655_s1605_typed_nadsoliton_metric_state_space_scale_quotient_pretheorem.py"
OUT = ROOT / "generated" / "p2655_s1605_typed_nadsoliton_metric_state_space_scale_quotient_pretheorem.json"
MD = ROOT / "generated" / "p2655_s1605_typed_nadsoliton_metric_state_space_scale_quotient_pretheorem.md"


class P2655TypedMetricScaleQuotientPretheoremTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        subprocess.run([sys.executable, str(SCRIPT)], cwd=ROOT, check=True)
        cls.payload = json.loads(OUT.read_text(encoding="utf-8"))

    def test_content_first_grep_and_upstream_consistency(self) -> None:
        audit = self.payload["semantic_rg_antiduplication_audit"]
        self.assertIn("content-first", audit["mode"])
        self.assertIn("typed_state_space_metric_content", audit["patterns"])
        self.assertIn("scale_quotient_content", audit["patterns"])
        self.assertIn("operator_identity_content", audit["patterns"])
        self.assertTrue(all(data["count"] > 0 for data in audit["patterns"].values()))
        upstream = self.payload["upstream_consistency"]
        self.assertTrue(upstream["p2653_typed_metric_uv_source_not_exported"])
        self.assertTrue(upstream["p2653_missing_atoms_include_uv_unit"])
        self.assertTrue(upstream["p2654_scale_invariant_selector_blocked"])
        self.assertTrue(upstream["p2654_raw_anchor_not_promoted"])

    def test_metric_skeleton_passes_but_does_not_select_unit(self) -> None:
        metric = self.payload["metric_axiom_audit"]
        self.assertTrue(metric["all_scales_metric_axioms_pass"])
        self.assertEqual(metric["normalized_geometry_rank_on_scale_orbit"], 1)
        self.assertLess(metric["max_normalized_difference"], 1e-12)
        self.assertFalse(metric["unit_selected_by_metric_axioms_alone"])
        self.assertTrue(self.payload["typed_state_space"]["no_sub_nadsoliton_information_layer"])

    def test_damping_covariance_is_verified_but_not_beta_source(self) -> None:
        covariance = self.payload["damping_covariance_audit"]
        self.assertTrue(covariance["covariance_verified"])
        self.assertLess(covariance["max_covariance_error"], 1e-12)
        self.assertFalse(covariance["beta_one_selected_by_covariance"])

    def test_source_matrix_blocks_current_candidates(self) -> None:
        matrix = {row["candidate"]: row for row in self.payload["source_candidate_matrix"]}
        self.assertFalse(matrix["typed_state_space_plus_metric_axioms_only"]["passes_as_typed_metric_uv_source_now"])
        self.assertFalse(matrix["typed_state_space_plus_metric_axioms_only"]["breaks_scale_orbit"])
        self.assertTrue(matrix["diameter_equals_one_normalization"]["uses_external_unit_anchor"])
        self.assertFalse(matrix["diameter_equals_one_normalization"]["passes_as_typed_metric_uv_source_now"])
        self.assertFalse(matrix["strict_denominator_covariance_identity"]["passes_as_typed_metric_uv_source_now"])
        self.assertFalse(matrix["typed_nadsoliton_dynamics_plus_conserved_uv_action_quantum"]["passes_as_typed_metric_uv_source_now"])

    def test_no_closure_and_docs_updated(self) -> None:
        decision = self.payload["closure_decision"]
        self.assertEqual(decision["passing_typed_metric_uv_source_candidates"], [])
        self.assertTrue(decision["typed_metric_skeleton_constructed"])
        self.assertTrue(decision["scale_quotient_rank_one"])
        self.assertTrue(decision["strict_denominator_covariance_verified"])
        self.assertFalse(decision["uv_unit_selected_now"])
        self.assertFalse(decision["beta_source_exported_now"])
        self.assertFalse(decision["role_bearing_ltotal_now"])
        self.assertFalse(decision["toe_closure_now"])
        self.assertTrue(all(value is False for value in self.payload["negative_export_flags"].values()))
        self.assertIn("P2655/S1605", MD.read_text(encoding="utf-8"))
        self.assertIn("P2655/S1605", (ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md").read_text(encoding="utf-8"))
        self.assertIn("P2655/S1605", (ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md").read_text(encoding="utf-8"))


if __name__ == "__main__":
    unittest.main()
