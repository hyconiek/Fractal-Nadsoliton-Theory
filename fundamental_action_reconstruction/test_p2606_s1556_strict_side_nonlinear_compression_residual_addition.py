from __future__ import annotations

import json
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT))

from p2606_s1556_strict_side_nonlinear_compression_residual_addition import MD, OUT, append_doc_sections, build_payload, write_markdown


class P2606StrictSideNonlinearCompressionResidualAdditionTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.payload = build_payload()
        OUT.write_text(json.dumps(cls.payload, indent=2, sort_keys=True), encoding="utf-8")
        write_markdown(cls.payload)
        append_doc_sections()
        cls.theorem = cls.payload["strict_side_nonlinear_compression_residual_addition"]["theorem_export"]
        cls.audit = cls.theorem["strict_side_nonlinear_residual_audit"]
        cls.remaining = cls.theorem["post_nonlinear_remaining_bridge_matrix"]

    def test_identity_and_inheritance(self) -> None:
        self.assertEqual(self.payload["packet_id"], "P2606")
        self.assertEqual(self.payload["stage_id"], "S1556")
        self.assertIn("NONLINEAR_COMPRESSION_RESIDUAL_COMPONENT", self.payload["status"])
        self.assertTrue(self.theorem["p2603_codimension_slope_inherited"])
        self.assertTrue(self.theorem["p2605_linear_slice_evidence_inherited"])

    def test_nonlinear_residual_component(self) -> None:
        self.assertTrue(self.theorem["nonlinear_compression_residual_addition_exported"])
        self.assertEqual(self.audit["eta_linear_slice"], 1.0)
        self.assertEqual(self.audit["eta_nonlinear_codimension"], 0.8)
        self.assertGreater(self.audit["max_abs_nonlinear_residual"], 1e-4)
        self.assertGreater(self.audit["nonzero_residual_count"], 0)
        self.assertTrue(self.audit["nonlinear_residual_not_absorbed_by_amplitude_rescaling"])
        self.assertGreater(self.audit["max_abs_residual_after_best_scalar_fit"], 1e-4)

    def test_remaining_bridge_matrix(self) -> None:
        self.assertEqual(self.remaining["truth_table_row_count"], 4)
        self.assertEqual(self.remaining["accepting_row_count"], 1)
        self.assertEqual(self.remaining["remaining_missing_gate_count_after_p2606"], 2)
        self.assertFalse(self.remaining["current_role_bearing_ltotal_accepts"])
        self.assertTrue(self.remaining["current_assignment_after_p2606"]["nonlinear_compression_residual_addition"])

    def test_scope_guards_and_docs(self) -> None:
        self.assertFalse(self.theorem["strict_side_residual_additions_exported"])
        self.assertFalse(self.theorem["phase_topological_selector_data_exported"])
        self.assertFalse(self.theorem["strict_damping_role_transfer_theorem_exported"])
        self.assertFalse(self.theorem["role_bearing_ltotal_exported"])
        self.assertFalse(self.theorem["qw2191_discharged_by_this_addition"])
        self.assertFalse(self.theorem["toe_closure_claimed"])
        self.assertTrue(all(self.payload["gatekeeper_checks"].values()))
        self.assertEqual(self.payload["rg_non_duplication_audit"]["tool"], "rg")
        self.assertIn("P2606/S1556", MD.read_text(encoding="utf-8"))
        self.assertIn("P2606/S1556", (ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md").read_text(encoding="utf-8"))
        self.assertIn("P2606/S1556", (ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md").read_text(encoding="utf-8"))


if __name__ == "__main__":
    unittest.main()
