from __future__ import annotations

import json
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT))

from p2604_s1554_strict_damping_post_source_bridge_readiness_matrix import MD, OUT, append_doc_sections, build_payload, write_markdown


class P2604StrictDampingPostSourceBridgeReadinessMatrixTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.payload = build_payload()
        OUT.write_text(json.dumps(cls.payload, indent=2, sort_keys=True), encoding="utf-8")
        write_markdown(cls.payload)
        append_doc_sections()
        cls.theorem = cls.payload["strict_damping_post_source_bridge_readiness_matrix"]["theorem_export"]

    def test_identity_and_p2603_inheritance(self) -> None:
        self.assertEqual(self.payload["packet_id"], "P2604")
        self.assertEqual(self.payload["stage_id"], "S1554")
        self.assertIn("POST_SOURCE_BRIDGE_READINESS_MATRIX", self.payload["status"])
        self.assertTrue(self.theorem["strict_damping_beta_eta_source_inherited_from_p2603"])
        self.assertTrue(self.payload["gatekeeper_checks"]["p2603_strict_damping_source_inherited"])

    def test_bridge_role_truth_table(self) -> None:
        self.assertEqual(self.theorem["bridge_readiness_truth_table_row_count"], 8)
        self.assertEqual(self.theorem["role_bearing_accepting_row_count"], 1)
        self.assertEqual(len(self.theorem["bridge_role_gates"]), 3)
        accepting = self.theorem["role_bearing_accepting_row"]
        self.assertTrue(accepting["role_bearing_ltotal_accepts"])
        self.assertTrue(all(accepting["assignment"].values()))
        self.assertFalse(self.theorem["current_role_bearing_ltotal_accepts"])
        self.assertEqual(len(self.theorem["concrete_missing_ingredients"]), 3)

    def test_scope_guards_and_docs(self) -> None:
        self.assertFalse(self.theorem["legacy_to_strict_completion_bridge_exported"])
        self.assertFalse(self.theorem["strict_side_residual_additions_exported"])
        self.assertFalse(self.theorem["strict_damping_role_transfer_theorem_exported"])
        self.assertFalse(self.theorem["role_bearing_ltotal_exported"])
        self.assertFalse(self.theorem["qw2191_discharged_by_this_matrix"])
        self.assertFalse(self.theorem["toe_closure_claimed"])
        self.assertTrue(all(self.payload["gatekeeper_checks"].values()))
        self.assertEqual(self.payload["rg_non_duplication_audit"]["tool"], "rg")
        self.assertIn("P2604/S1554", MD.read_text(encoding="utf-8"))
        self.assertIn("P2604/S1554", (ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md").read_text(encoding="utf-8"))
        self.assertIn("P2604/S1554", (ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md").read_text(encoding="utf-8"))


if __name__ == "__main__":
    unittest.main()
