from __future__ import annotations

import json
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT))

from p2497_s1447_phase_normalized_curvature_inflection_window_interval_isolation_certificate import (
    MD,
    OUT,
    build_payload,
    write_markdown,
)


class P2497PhaseNormalizedCurvatureInflectionWindowIntervalIsolationCertificateTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.payload = build_payload()
        OUT.write_text(json.dumps(cls.payload, indent=2, sort_keys=True), encoding="utf-8")
        write_markdown(cls.payload)
        cls.theorem = cls.payload["phase_normalized_curvature_inflection_window_interval_isolation_certificate"]["theorem_export"]
        cls.cert = cls.theorem["inflection_window_interval_isolation_certificate"]

    def test_identity_and_status(self) -> None:
        self.assertEqual(self.payload["packet_id"], "P2497")
        self.assertEqual(self.payload["stage_id"], "S1447")
        self.assertIn("INFLECTION_WINDOW_INTERVAL_ISOLATION", self.payload["status"])
        self.assertIn("NO_FORMAL_DIRECTED_BACKEND", self.payload["status"])
        self.assertIn("NO_GLOBAL_UNIQUENESS", self.payload["status"])
        self.assertIn("NO_SOURCE_EXPORT", self.payload["status"])
        self.assertIn("NO_BRIDGE_ATOM", self.payload["status"])
        self.assertIn("NO_QW2191", self.payload["status"])
        self.assertIn("NO_ROLE_TRANSFER", self.payload["status"])
        self.assertIn("NO_TOE", self.payload["status"])

    def test_root_window_and_endpoint_signs(self) -> None:
        self.assertEqual(self.theorem["interval_backend"], "mpmath.iv")
        self.assertEqual(self.theorem["interval_dps"], 90)
        self.assertEqual(self.theorem["point_solve_dps"], 100)
        self.assertEqual(self.theorem["inflection_root_window_d"], ["0.3495", "0.3498"])
        self.assertLessEqual(float(self.theorem["inflection_root_window_width"]), 0.0003)
        self.assertTrue(self.theorem["point_root_inside_window"])
        self.assertTrue(self.theorem["endpoint_intervals_have_opposite_strict_signs"])
        self.assertEqual(self.cert["window_left_endpoint_interval"]["curvature_interval_sign"], 1)
        self.assertEqual(self.cert["window_right_endpoint_interval"]["curvature_interval_sign"], -1)

    def test_adaptive_slab_exclusions(self) -> None:
        left = self.cert["left_positive_slab_exclusion"]
        right = self.cert["right_negative_slab_exclusion"]
        self.assertTrue(left["all_cells_exclude_zero_with_expected_sign"])
        self.assertTrue(right["all_cells_exclude_zero_with_expected_sign"])
        self.assertEqual(left["unresolved_cell_count"], 0)
        self.assertEqual(right["unresolved_cell_count"], 0)
        self.assertGreater(self.theorem["left_positive_slab_cell_count"], 0)
        self.assertGreater(self.theorem["right_negative_slab_cell_count"], 0)
        self.assertTrue(all(row["curvature_interval_sign"] == 1 for row in left["accepted_cells"]))
        self.assertTrue(all(row["curvature_interval_sign"] == -1 for row in right["accepted_cells"]))
        self.assertTrue(self.theorem["outside_window_interval_zero_exclusion_certified"])
        self.assertTrue(self.theorem["single_unresolved_inflection_window_on_audited_domain"])

    def test_inherited_branch_guard_and_gatekeepers(self) -> None:
        self.assertTrue(self.theorem["p2496_inverse_branch_inherited"])
        self.assertTrue(all(self.payload["gatekeeper_checks"].values()))

    def test_negative_controls(self) -> None:
        self.assertFalse(self.theorem["formal_directed_rounding_backend_exported"])
        self.assertFalse(self.theorem["global_inflection_uniqueness_theorem_exported"])
        self.assertFalse(self.theorem["analytic_monotonicity_theorem_exported"])
        self.assertFalse(self.theorem["curvature_dynamic_source_exported"])
        self.assertFalse(self.theorem["legacy_to_strict_bridge_atom_exported"])
        self.assertFalse(self.theorem["strict_compression_dynamic_source_exported"])
        self.assertFalse(self.theorem["selector_source_theorem_exported"])
        self.assertFalse(self.theorem["qw2191_discharged_by_this_certificate"])
        self.assertFalse(self.theorem["role_transfer_licensed_by_this_certificate"])
        self.assertFalse(self.theorem["toe_closure_exported"])

    def test_rg_audit_and_documentation_markers(self) -> None:
        audit = self.payload["rg_non_duplication_audit"]
        self.assertEqual(audit["tool"], "rg")
        self.assertIn("new_packet", audit["patterns"])
        self.assertGreaterEqual(audit["patterns"]["new_packet"]["count"], 0)
        self.assertIn("P2497/S1447", MD.read_text(encoding="utf-8"))
        self.assertIn("P2497/S1447", (ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md").read_text(encoding="utf-8"))
        self.assertIn("P2497/S1447", (ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md").read_text(encoding="utf-8"))


if __name__ == "__main__":
    unittest.main()
