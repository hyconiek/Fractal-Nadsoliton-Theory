from __future__ import annotations

import json
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT))

from p2554_s1504_strict_damping_local_exhaustion_bridge_reorientation_certificate import GATES, MD, OUT, append_doc_sections, build_payload, write_markdown


class P2554LocalExhaustionBridgeReorientationTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.payload = build_payload()
        OUT.write_text(json.dumps(cls.payload, indent=2, sort_keys=True), encoding="utf-8")
        write_markdown(cls.payload)
        append_doc_sections()
        cls.theorem = cls.payload["strict_damping_local_exhaustion_bridge_reorientation_certificate"]["theorem_export"]

    def test_identity_and_precursors(self) -> None:
        self.assertEqual(self.payload["packet_id"], "P2554")
        self.assertEqual(self.payload["stage_id"], "S1504")
        self.assertIn("LOCAL_EXHAUSTION_BRIDGE_REORIENTATION", self.payload["status"])
        self.assertTrue(self.theorem["p2539_toe_potential_inherited"])
        self.assertTrue(self.theorem["p2547_post_identity_residual_trikey_inherited"])
        self.assertTrue(self.theorem["p2553_anchor_constant_obligation_inherited"])

    def test_bridge_gate_truth_table(self) -> None:
        self.assertEqual(self.theorem["bridge_gate_vector"], GATES)
        self.assertEqual(self.theorem["bridge_gate_truth_table_row_count"], 64)
        self.assertEqual(self.theorem["bridge_gate_accepting_row_count"], 1)
        self.assertTrue(all(self.theorem["bridge_gate_accepting_row"]["assignment"].values()))
        self.assertEqual(self.theorem["single_gate_sensitivity_row_count"], len(GATES))
        self.assertTrue(self.theorem["all_single_gate_omissions_reject"])
        self.assertFalse(self.theorem["local_strict_damping_source_alone_accepts_toe_or_role_transfer"])
        self.assertTrue(self.theorem["local_route_exhaustion_exported_as_reorientation_not_closure"])

    def test_negative_controls_and_recommendation(self) -> None:
        self.assertFalse(self.theorem["strict_damping_local_source_obligation_discharge_exported"])
        self.assertFalse(self.theorem["legacy_to_strict_completion_bridge_exported"])
        self.assertFalse(self.theorem["role_transfer_theorem_exported"])
        self.assertFalse(self.theorem["qw2191_selector_discharge_exported"])
        self.assertFalse(self.theorem["role_bearing_ltotal_exported"])
        self.assertFalse(self.theorem["toe_closure_claimed"])
        self.assertFalse(self.theorem["legacy_role_transfer_claimed"])
        self.assertIn("legacy->strict completion/source bridge audit", self.theorem["recommended_next_honest_step"])
        self.assertTrue(all(self.payload["gatekeeper_checks"].values()))

    def test_rg_and_docs(self) -> None:
        self.assertEqual(self.payload["rg_non_duplication_audit"]["tool"], "rg")
        self.assertIn("P2554/S1504", MD.read_text(encoding="utf-8"))
        self.assertIn("P2554/S1504", (ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md").read_text(encoding="utf-8"))
        self.assertIn("P2554/S1504", (ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md").read_text(encoding="utf-8"))


if __name__ == "__main__":
    unittest.main()
