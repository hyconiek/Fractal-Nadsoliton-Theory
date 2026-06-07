from __future__ import annotations

import json
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT))

from p2560_s1510_legacy_to_strict_damping_constant_log_source_current_premise_obstruction_certificate import MD, OUT, append_doc_sections, build_payload, write_markdown


class P2560ConstantLogSourceCurrentPremiseObstructionTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.payload = build_payload()
        OUT.write_text(json.dumps(cls.payload, indent=2, sort_keys=True), encoding="utf-8")
        write_markdown(cls.payload)
        append_doc_sections()
        cls.theorem = cls.payload["legacy_to_strict_damping_constant_log_source_current_premise_obstruction_certificate"]["theorem_export"]

    def test_countermodels_block_current_premise_entailment(self) -> None:
        self.assertEqual(self.payload["packet_id"], "P2560")
        self.assertEqual(self.payload["stage_id"], "S1510")
        self.assertTrue(self.theorem["p2559_conditional_selector_inherited"])
        self.assertEqual(self.theorem["countermodel_count"], 88)
        self.assertEqual(self.theorem["countermodel_q_values"], [-2.0, -1.0, 1.0, 2.0])
        self.assertTrue(self.theorem["all_countermodels_pass_endpoint_premises"])
        self.assertTrue(self.theorem["all_countermodels_have_positive_midpoint_denominator"])
        self.assertTrue(self.theorem["all_countermodels_violate_constant_log_source"])
        self.assertTrue(self.theorem["current_endpoint_premises_do_not_entail_constant_log_source"])
        self.assertTrue(self.theorem["p2559_conditional_selector_cannot_be_promoted_to_source_by_current_premises"])

    def test_negative_controls_and_docs(self) -> None:
        self.assertFalse(self.theorem["constant_log_source_law_exported"])
        self.assertFalse(self.theorem["geometric_homotopy_source_exported"])
        self.assertFalse(self.theorem["unique_damping_homotopy_source_exported"])
        self.assertFalse(self.theorem["damping_compression_bridge_component_ready"])
        self.assertFalse(self.theorem["full_bridge_theorem_exported"])
        self.assertFalse(self.theorem["role_transfer_theorem_exported"])
        self.assertFalse(self.theorem["toe_closure_claimed"])
        self.assertIn("derive a strict dynamic source", self.theorem["recommended_next_honest_step"])
        self.assertTrue(all(self.payload["gatekeeper_checks"].values()))
        self.assertEqual(self.payload["rg_non_duplication_audit"]["tool"], "rg")
        self.assertIn("P2560/S1510", MD.read_text(encoding="utf-8"))
        self.assertIn("P2560/S1510", (ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md").read_text(encoding="utf-8"))
        self.assertIn("P2560/S1510", (ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md").read_text(encoding="utf-8"))


if __name__ == "__main__":
    unittest.main()
