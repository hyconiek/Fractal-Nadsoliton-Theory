from __future__ import annotations

import json
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT))

from p2611_s1561_ltotal_role_semantics_acceptance_predicate import ACCEPTANCE_GATES, MD, OUT, append_doc_sections, build_payload, write_markdown


class P2611LtotalRoleSemanticsAcceptancePredicateTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.payload = build_payload()
        OUT.write_text(json.dumps(cls.payload, indent=2, sort_keys=True), encoding="utf-8")
        write_markdown(cls.payload)
        append_doc_sections()
        cls.theorem = cls.payload["ltotal_role_semantics_acceptance_predicate"]["theorem_export"]
        cls.table = cls.theorem["acceptance_truth_table"]

    def test_identity_and_semantics(self) -> None:
        self.assertEqual(self.payload["packet_id"], "P2611")
        self.assertEqual(self.payload["stage_id"], "S1561")
        self.assertIn("ROLE_SEMANTICS_DEFINED", self.payload["status"])
        self.assertTrue(self.theorem["role_semantics_definition_exported"])
        self.assertEqual(self.theorem["acceptance_gates"], ACCEPTANCE_GATES)
        self.assertEqual(len(self.theorem["role_semantics_axioms"]), 6)

    def test_acceptance_truth_table_rejects_current_assignment(self) -> None:
        self.assertEqual(self.table["truth_table_row_count"], 16)
        self.assertEqual(self.table["accepting_row_count"], 1)
        self.assertFalse(self.table["current_assignment_accepts"])
        self.assertTrue(self.table["current_assignment"]["formal_role_semantics_defined"])
        self.assertFalse(self.table["current_assignment"]["strict_damping_source_revalidated_after_p2610"])
        self.assertFalse(self.table["current_assignment"]["bridge_completion_revalidated_after_p2610"])
        self.assertFalse(self.table["current_assignment"]["role_transfer_revalidated_after_p2610"])

    def test_p2610_quarantine_and_no_reexports(self) -> None:
        self.assertTrue(self.theorem["p2610_quarantine_inherited"]["p2607_bridge_quarantined"])
        self.assertTrue(self.theorem["p2610_quarantine_inherited"]["p2608_role_transfer_quarantined"])
        self.assertFalse(self.theorem["strict_damping_beta_eta_source_reexported"])
        self.assertFalse(self.theorem["bridge_completion_reexported"])
        self.assertFalse(self.theorem["role_transfer_theorem_reexported"])
        self.assertFalse(self.theorem["role_bearing_ltotal_accepted"])
        self.assertFalse(self.theorem["qw2191_discharged_by_this_packet"])
        self.assertFalse(self.theorem["toe_closure_claimed"])
        self.assertFalse(self.theorem["apd_source_exported"])

    def test_gatekeepers_and_docs(self) -> None:
        self.assertTrue(all(self.payload["gatekeeper_checks"].values()))
        self.assertEqual(self.payload["rg_non_duplication_audit"]["tool"], "rg")
        self.assertIn("P2611/S1561", MD.read_text(encoding="utf-8"))
        self.assertIn("P2611/S1561", (ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md").read_text(encoding="utf-8"))
        self.assertIn("P2611/S1561", (ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md").read_text(encoding="utf-8"))


if __name__ == "__main__":
    unittest.main()
