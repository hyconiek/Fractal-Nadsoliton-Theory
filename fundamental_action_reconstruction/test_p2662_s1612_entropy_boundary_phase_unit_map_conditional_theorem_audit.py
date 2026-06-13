from __future__ import annotations

import json
import subprocess
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent
SCRIPT = ROOT / "p2662_s1612_entropy_boundary_phase_unit_map_conditional_theorem_audit.py"
OUT = ROOT / "generated" / "p2662_s1612_entropy_boundary_phase_unit_map_conditional_theorem_audit.json"
MD = ROOT / "generated" / "p2662_s1612_entropy_boundary_phase_unit_map_conditional_theorem_audit.md"


class P2662EntropyBoundaryPhaseUnitMapConditionalTheoremAuditTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        subprocess.run([sys.executable, str(SCRIPT)], cwd=ROOT, check=True)
        cls.payload = json.loads(OUT.read_text(encoding="utf-8"))

    def test_content_first_and_upstream_consistency(self) -> None:
        audit = self.payload["semantic_rg_antiduplication_audit"]
        self.assertIn("content-first", audit["mode"])
        for key in (
            "entropy_measure_content",
            "boundary_phase_content",
            "bit_action_length_content",
            "nonclosure_guard_content",
        ):
            self.assertIn(key, audit["patterns"])
            self.assertGreater(audit["patterns"][key]["count"], 0)
        upstream = self.payload["upstream_consistency"]
        self.assertTrue(upstream["p2660_all_topological_candidates_need_unit_map"])
        self.assertTrue(upstream["p2660_no_beta_source"])
        self.assertTrue(upstream["p2661_fractal_log_shift_verified"])
        self.assertTrue(upstream["p2661_selection_not_intrinsic_without_reference_measure"])
        self.assertTrue(upstream["p2661_qw2191_not_discharged"])

    def test_conditional_theorem_selects_unique_covariant_physical_scale(self) -> None:
        conditional = self.payload["conditional_unit_map_theorem_audit"]
        self.assertTrue(conditional["conditional_unique_positive_scale_for_each_integer_source"])
        self.assertTrue(conditional["physical_scale_invariant_under_base_coordinate_rescaling"])
        self.assertLess(conditional["max_entropy_residual"], 1e-11)
        self.assertFalse(conditional["premise_intrinsic_entropy_measure_exported"])
        self.assertFalse(conditional["premise_boundary_phase_bit_target_exported"])
        self.assertFalse(conditional["premise_bit_to_action_or_length_map_exported"])
        self.assertFalse(conditional["unconditional_uv_unit_selected_now"])

    def test_premise_gap_ledger_keeps_theorem_conditional(self) -> None:
        gaps = {gap["premise"]: gap for gap in self.payload["premise_gap_audit"]}
        for premise in (
            "intrinsic_pre_normalization_entropy_measure",
            "boundary_phase_integer_to_entropy_target",
            "bit_to_action_or_bit_to_length_unit_map",
            "selector_branch_orientation_law",
        ):
            self.assertIn(premise, gaps)
            self.assertFalse(gaps[premise]["currently_exported"])

    def test_source_matrix_blocks_unconditional_false_pass(self) -> None:
        matrix = {row["candidate"]: row for row in self.payload["source_candidate_matrix"]}
        self.assertTrue(matrix["conditional_entropy_boundary_phase_unit_map"]["selects_unique_scale_conditionally"])
        self.assertFalse(matrix["conditional_entropy_boundary_phase_unit_map"]["all_premises_exported"])
        self.assertFalse(matrix["conditional_entropy_boundary_phase_unit_map"]["passes_as_uv_unit_source_now"])
        self.assertFalse(matrix["raw_cocycle_integer_as_entropy_target"]["passes_as_uv_unit_source_now"])
        self.assertFalse(matrix["bit_to_action_or_length_conversion"]["passes_as_uv_unit_source_now"])
        self.assertFalse(matrix["entropy_selector_branch_for_qw2191"]["passes_as_uv_unit_source_now"])

    def test_no_closure_and_docs_updated(self) -> None:
        decision = self.payload["closure_decision"]
        self.assertEqual(decision["passing_uv_unit_source_candidates"], [])
        self.assertTrue(decision["conditional_unique_scale_verified"])
        self.assertTrue(decision["physical_scale_covariance_verified"])
        self.assertFalse(decision["unconditional_uv_unit_selected_now"])
        self.assertFalse(decision["beta_source_exported_now"])
        self.assertFalse(decision["qw2191_discharged_now"])
        self.assertFalse(decision["role_bearing_ltotal_now"])
        self.assertFalse(decision["toe_closure_now"])
        self.assertTrue(all(value is False for value in self.payload["negative_export_flags"].values()))
        self.assertIn("P2662/S1612", MD.read_text(encoding="utf-8"))
        self.assertIn("P2662/S1612", (ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md").read_text(encoding="utf-8"))
        self.assertIn("P2662/S1612", (ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md").read_text(encoding="utf-8"))


if __name__ == "__main__":
    unittest.main()
