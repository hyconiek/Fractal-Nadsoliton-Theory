from __future__ import annotations

import json
import subprocess
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent
SCRIPT = ROOT / "p2661_s1611_shannon_entropy_scale_anomaly_uv_anchor_audit.py"
OUT = ROOT / "generated" / "p2661_s1611_shannon_entropy_scale_anomaly_uv_anchor_audit.json"
MD = ROOT / "generated" / "p2661_s1611_shannon_entropy_scale_anomaly_uv_anchor_audit.md"


class P2661ShannonEntropyScaleAnomalyUvAnchorAuditTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        subprocess.run([sys.executable, str(SCRIPT)], cwd=ROOT, check=True)
        cls.payload = json.loads(OUT.read_text(encoding="utf-8"))

    def test_content_first_and_upstream_consistency(self) -> None:
        audit = self.payload["semantic_rg_antiduplication_audit"]
        self.assertIn("content-first", audit["mode"])
        for key in (
            "shannon_entropy_content",
            "scale_anomaly_content",
            "uv_anchor_content",
            "nonclosure_guard_content",
        ):
            self.assertIn(key, audit["patterns"])
            self.assertGreater(audit["patterns"][key]["count"], 0)
        upstream = self.payload["upstream_consistency"]
        self.assertTrue(upstream["p2658_homogeneous_no_go_verified"])
        self.assertTrue(upstream["p2659_no_intrinsic_anomaly_source"])
        self.assertTrue(upstream["p2660_requires_dimensionful_unit_map"])
        self.assertTrue(upstream["p2660_no_beta_source"])

    def test_normalized_shannon_entropy_cancels_global_scale(self) -> None:
        normalized = self.payload["normalized_shannon_audit"]
        self.assertTrue(normalized["normalized_entropy_is_scale_invariant"])
        self.assertLess(normalized["max_entropy_drift"], 1e-11)
        self.assertLess(normalized["max_distribution_drift"], 1e-11)
        self.assertFalse(normalized["normalized_entropy_selects_unique_scale"])

    def test_fractal_entropy_has_log_shift_but_needs_reference(self) -> None:
        fractal = self.payload["fractal_differential_entropy_audit"]
        self.assertTrue(fractal["additive_log_scale_shift_verified"])
        self.assertLess(fractal["max_shift_error"], 1e-11)
        self.assertTrue(fractal["declared_bit_levels_can_select_representatives"])
        self.assertFalse(fractal["selection_is_intrinsic_without_reference_measure"])
        self.assertTrue(fractal["bit_is_dimensionless_not_length_unit"])

    def test_entropy_arrow_does_not_discharge_qw2191(self) -> None:
        arrow = self.payload["entropy_arrow_audit"]
        self.assertTrue(arrow["entropy_increases_with_scale_in_model"])
        self.assertFalse(arrow["exports_time_arrow_source"])
        self.assertFalse(arrow["discharges_qw_2191"])

    def test_source_matrix_blocks_entropy_false_passes(self) -> None:
        matrix = {row["candidate"]: row for row in self.payload["source_candidate_matrix"]}
        self.assertFalse(matrix["normalized_finite_shannon_entropy_of_distance_weights"]["has_log_scale_anomaly"])
        self.assertFalse(matrix["normalized_finite_shannon_entropy_of_distance_weights"]["passes_as_uv_unit_source_now"])
        self.assertTrue(matrix["fractal_differential_entropy_Df_log_a_shift"]["has_log_scale_anomaly"])
        self.assertTrue(matrix["fractal_differential_entropy_Df_log_a_shift"]["uses_external_reference_or_unit"])
        self.assertFalse(matrix["fractal_differential_entropy_Df_log_a_shift"]["passes_as_uv_unit_source_now"])
        self.assertFalse(matrix["one_bit_entropy_quantum_as_length_unit"]["passes_as_uv_unit_source_now"])
        self.assertFalse(matrix["entropy_arrow_as_qw2191_selector"]["passes_as_uv_unit_source_now"])
        self.assertFalse(matrix["intrinsic_nadsoliton_entropy_measure_and_unit_map_theorem"]["covered_by_audit"])

    def test_no_closure_and_docs_updated(self) -> None:
        decision = self.payload["closure_decision"]
        self.assertEqual(decision["passing_uv_unit_source_candidates"], [])
        self.assertTrue(decision["normalized_entropy_is_scale_invariant"])
        self.assertTrue(decision["fractal_log_shift_verified"])
        self.assertFalse(decision["selection_is_intrinsic_without_reference_measure"])
        self.assertFalse(decision["entropy_arrow_discharges_qw_2191"])
        self.assertFalse(decision["uv_unit_selected_now"])
        self.assertFalse(decision["beta_source_exported_now"])
        self.assertFalse(decision["role_bearing_ltotal_now"])
        self.assertFalse(decision["toe_closure_now"])
        self.assertTrue(all(value is False for value in self.payload["negative_export_flags"].values()))
        self.assertIn("P2661/S1611", MD.read_text(encoding="utf-8"))
        self.assertIn("P2661/S1611", (ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md").read_text(encoding="utf-8"))
        self.assertIn("P2661/S1611", (ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md").read_text(encoding="utf-8"))


if __name__ == "__main__":
    unittest.main()
