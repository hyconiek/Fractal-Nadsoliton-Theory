from __future__ import annotations

import json
import subprocess
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent
SCRIPT = ROOT / "p2660_s1610_boundary_cocycle_anomaly_coefficient_dimension_audit.py"
OUT = ROOT / "generated" / "p2660_s1610_boundary_cocycle_anomaly_coefficient_dimension_audit.json"
MD = ROOT / "generated" / "p2660_s1610_boundary_cocycle_anomaly_coefficient_dimension_audit.md"


class P2660BoundaryCocycleAnomalyCoefficientDimensionAuditTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        subprocess.run([sys.executable, str(SCRIPT)], cwd=ROOT, check=True)
        cls.payload = json.loads(OUT.read_text(encoding="utf-8"))

    def test_content_first_and_upstream_consistency(self) -> None:
        audit = self.payload["semantic_rg_antiduplication_audit"]
        self.assertIn("content-first", audit["mode"])
        for key in (
            "boundary_cocycle_content",
            "anomaly_coefficient_content",
            "dimension_unit_content",
            "nonclosure_guard_content",
        ):
            self.assertIn(key, audit["patterns"])
            self.assertGreater(audit["patterns"][key]["count"], 0)
        upstream = self.payload["upstream_consistency"]
        self.assertTrue(upstream["p2658_homogeneous_no_go_verified"])
        self.assertTrue(upstream["p2659_no_intrinsic_anomaly_source"])
        self.assertTrue(upstream["p2659_next_step_requests_boundary_cocycle_phase_law"])

    def test_boundary_cocycle_invariants_are_scale_invariant_but_not_units(self) -> None:
        invariants = self.payload["boundary_cocycle_invariants"]
        self.assertEqual(invariants["graph_beta_1"], 2)
        self.assertEqual(invariants["beta_1_after_triangle_fill"], 1)
        self.assertEqual(invariants["euler_characteristic"], 0)
        self.assertTrue(invariants["topological_numbers_are_scale_invariant"])
        dimension = self.payload["dimension_typing_audit"]
        self.assertTrue(dimension["all_topological_candidates_need_unit_map"])
        self.assertFalse(dimension["dimensionful_anomaly_coefficient_derived_now"])
        self.assertTrue(all(not row["can_be_added_to_trace_action_without_unit_map"] for row in dimension["rows"]))

    def test_topological_phase_audit_still_has_clock_freedom(self) -> None:
        phase = self.payload["topological_anomaly_phase_audit"]
        self.assertTrue(phase["all_integer_phase_conditions_satisfied"])
        self.assertLess(phase["max_phase_error"], 1e-11)
        self.assertTrue(phase["declared_quantum_selectors_are_external"])
        self.assertFalse(phase["uv_unit_selected_by_topological_anomaly_now"])

    def test_source_matrix_blocks_topological_false_pass(self) -> None:
        matrix = {row["candidate"]: row for row in self.payload["source_candidate_matrix"]}
        self.assertFalse(matrix["raw_boundary_cocycle_integer_as_lambda"]["dimensionally_admissible_without_unit_map"])
        self.assertFalse(matrix["raw_boundary_cocycle_integer_as_lambda"]["passes_as_uv_unit_source_now"])
        self.assertFalse(matrix["topological_integer_plus_integer_phase"]["passes_as_uv_unit_source_now"])
        self.assertTrue(matrix["declared_absolute_action_quantum_after_topology"]["uses_external_clock_or_scale_anchor"])
        self.assertFalse(matrix["declared_absolute_action_quantum_after_topology"]["passes_as_uv_unit_source_now"])
        self.assertFalse(matrix["derived_unit_map_from_nadsoliton_boundary_phase_law"]["covered_by_audit"])

    def test_no_closure_and_docs_updated(self) -> None:
        decision = self.payload["closure_decision"]
        self.assertEqual(decision["passing_uv_unit_source_candidates"], [])
        self.assertTrue(decision["all_topological_candidates_need_unit_map"])
        self.assertFalse(decision["uv_unit_selected_now"])
        self.assertFalse(decision["beta_source_exported_now"])
        self.assertFalse(decision["role_bearing_ltotal_now"])
        self.assertFalse(decision["toe_closure_now"])
        self.assertTrue(all(value is False for value in self.payload["negative_export_flags"].values()))
        self.assertIn("P2660/S1610", MD.read_text(encoding="utf-8"))
        self.assertIn("P2660/S1610", (ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md").read_text(encoding="utf-8"))
        self.assertIn("P2660/S1610", (ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md").read_text(encoding="utf-8"))


if __name__ == "__main__":
    unittest.main()
