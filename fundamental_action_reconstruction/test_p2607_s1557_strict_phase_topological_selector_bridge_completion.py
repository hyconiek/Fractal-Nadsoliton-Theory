from __future__ import annotations

import json
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT))

from p2607_s1557_strict_phase_topological_selector_bridge_completion import MD, OUT, append_doc_sections, build_payload, write_markdown


class P2607StrictPhaseTopologicalSelectorBridgeCompletionTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.payload = build_payload()
        OUT.write_text(json.dumps(cls.payload, indent=2, sort_keys=True), encoding="utf-8")
        write_markdown(cls.payload)
        append_doc_sections()
        cls.theorem = cls.payload["strict_phase_topological_selector_bridge_completion"]["theorem_export"]
        cls.audit = cls.theorem["phase_topological_selector_audit"]
        cls.role = cls.theorem["post_bridge_role_transfer_matrix"]

    def test_identity_and_inheritance(self) -> None:
        self.assertEqual(self.payload["packet_id"], "P2607")
        self.assertEqual(self.payload["stage_id"], "S1557")
        self.assertIn("PHASE_TOPOLOGICAL_SELECTOR_DATA", self.payload["status"])
        self.assertTrue(self.theorem["p2605_linear_slice_evidence_inherited"])
        self.assertTrue(self.theorem["p2606_nonlinear_residual_inherited"])

    def test_gf2_phase_selector_uniqueness(self) -> None:
        self.assertEqual(self.audit["rank_over_gf2"], 11)
        self.assertEqual(self.audit["nullity_over_gf2"], 0)
        self.assertTrue(self.audit["rank_11_uniqueness"])
        self.assertTrue(self.audit["all_single_bit_toggles_rejected"])
        self.assertEqual(len(self.audit["single_bit_toggle_witnesses"]), 11)
        self.assertTrue(self.theorem["phase_topological_selector_data_exported"])

    def test_bridge_completion_but_role_transfer_blocked(self) -> None:
        self.assertTrue(self.theorem["strict_side_residual_additions_exported"])
        self.assertTrue(self.theorem["legacy_to_strict_completion_bridge_exported"])
        self.assertEqual(self.role["truth_table_row_count"], 2)
        self.assertEqual(self.role["accepting_row_count"], 1)
        self.assertEqual(self.role["remaining_missing_gate_count_after_p2607"], 1)
        self.assertFalse(self.role["current_role_bearing_ltotal_accepts"])
        self.assertFalse(self.theorem["strict_damping_role_transfer_theorem_exported"])

    def test_scope_guards_and_docs(self) -> None:
        self.assertFalse(self.theorem["legacy_physical_role_transfer_exported"])
        self.assertFalse(self.theorem["role_bearing_ltotal_exported"])
        self.assertFalse(self.theorem["qw2191_discharged_by_this_completion"])
        self.assertFalse(self.theorem["toe_closure_claimed"])
        self.assertTrue(all(self.payload["gatekeeper_checks"].values()))
        self.assertEqual(self.payload["rg_non_duplication_audit"]["tool"], "rg")
        self.assertIn("P2607/S1557", MD.read_text(encoding="utf-8"))
        self.assertIn("P2607/S1557", (ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md").read_text(encoding="utf-8"))
        self.assertIn("P2607/S1557", (ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md").read_text(encoding="utf-8"))


if __name__ == "__main__":
    unittest.main()
