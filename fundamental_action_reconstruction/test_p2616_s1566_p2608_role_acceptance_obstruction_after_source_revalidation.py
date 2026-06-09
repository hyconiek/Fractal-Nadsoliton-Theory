from __future__ import annotations

import json
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT))

from p2616_s1566_p2608_role_acceptance_obstruction_after_source_revalidation import MD, OUT, append_doc_sections, build_payload, write_markdown


class P2616P2608RoleAcceptanceObstructionAfterSourceRevalidationTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.payload = build_payload()
        OUT.write_text(json.dumps(cls.payload, indent=2, sort_keys=True), encoding="utf-8")
        write_markdown(cls.payload)
        append_doc_sections()
        cls.theorem = cls.payload["p2608_role_acceptance_obstruction_after_source_revalidation"]["theorem_export"]

    def test_identity_and_truth_table(self) -> None:
        self.assertEqual(self.payload["packet_id"], "P2616")
        self.assertEqual(self.payload["stage_id"], "S1566")
        self.assertIn("ROLE_ACCEPTANCE_REJECTED", self.payload["status"])
        self.assertEqual(self.theorem["truth_table"]["row_count"], 16)
        self.assertEqual(self.theorem["truth_table"]["accepting_row_count"], 1)

    def test_source_repaired_but_bridge_blocks_acceptance(self) -> None:
        assignment = self.theorem["current_assignment_after_p2615"]
        self.assertTrue(assignment["formal_role_semantics_defined"])
        self.assertTrue(assignment["strict_damping_source_revalidated"])
        self.assertFalse(assignment["legacy_to_strict_bridge_revalidated"])
        self.assertFalse(assignment["role_transfer_revalidated"])
        self.assertFalse(self.theorem["current_assignment_accepts_role_bearing_ltotal"])
        self.assertEqual(self.theorem["current_missing_gates_after_p2615"], ["legacy_to_strict_bridge_revalidated", "role_transfer_revalidated"])

    def test_formal_obstruction_mentions_sources_and_bridge(self) -> None:
        proof_text = " ".join(self.theorem["formal_obstruction_proof"]["proof_steps"])
        self.assertIn("P2613", proof_text)
        self.assertIn("P2614", proof_text)
        self.assertIn("P2612", proof_text)
        self.assertIn("P2615", proof_text)
        self.assertTrue(self.theorem["p2608_quarantine_retained_by_p2616"])

    def test_scope_guards(self) -> None:
        self.assertFalse(self.theorem["p2608_role_bearing_ltotal_reenabled"])
        self.assertFalse(self.theorem["role_transfer_theorem_exported"])
        self.assertFalse(self.theorem["p2607_bridge_completion_revalidated"])
        self.assertFalse(self.theorem["legacy_physical_role_transfer_exported"])
        self.assertFalse(self.theorem["qw2191_discharged_by_this_packet"])
        self.assertFalse(self.theorem["toe_closure_claimed"])
        self.assertFalse(self.theorem["apd_source_exported"])

    def test_gatekeepers_and_docs(self) -> None:
        self.assertTrue(all(self.payload["gatekeeper_checks"].values()))
        self.assertEqual(self.payload["rg_non_duplication_audit"]["tool"], "rg")
        self.assertIn("P2616/S1566", MD.read_text(encoding="utf-8"))
        self.assertIn("P2616/S1566", (ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md").read_text(encoding="utf-8"))
        self.assertIn("P2616/S1566", (ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md").read_text(encoding="utf-8"))


if __name__ == "__main__":
    unittest.main()
