from __future__ import annotations

import json
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT))

from p2608_s1558_strict_damping_role_transfer_to_ltotal_theorem import MD, OUT, append_doc_sections, build_payload, write_markdown


class P2608StrictDampingRoleTransferToLtotalTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.payload = build_payload()
        OUT.write_text(json.dumps(cls.payload, indent=2, sort_keys=True), encoding="utf-8")
        write_markdown(cls.payload)
        append_doc_sections()
        cls.theorem = cls.payload["strict_damping_role_transfer_to_ltotal_theorem"]["theorem_export"]
        cls.table = cls.theorem["role_transfer_truth_table"]

    def test_identity_and_inheritance(self) -> None:
        self.assertEqual(self.payload["packet_id"], "P2608")
        self.assertEqual(self.payload["stage_id"], "S1558")
        self.assertIn("STRICT_DAMPING_ROLE_TRANSFER", self.payload["status"])
        self.assertTrue(self.theorem["strict_damping_beta_eta_source_inherited"])
        self.assertTrue(self.theorem["legacy_to_strict_completion_bridge_inherited"])

    def test_role_transfer_exports_scoped_ltotal(self) -> None:
        self.assertTrue(self.theorem["strict_damping_role_transfer_theorem_exported"])
        self.assertTrue(self.theorem["role_bearing_ltotal_exported"])
        self.assertEqual(self.table["truth_table_row_count"], 8)
        self.assertEqual(self.table["accepting_row_count"], 1)
        self.assertTrue(self.table["current_assignment_accepts"])
        self.assertEqual(self.theorem["role_bearing_ltotal_term"]["term_name"], "strict_damping_compression_beta_eta_term")

    def test_legacy_roles_remain_blocked(self) -> None:
        self.assertFalse(self.theorem["legacy_physical_role_transfer_exported"])
        self.assertTrue(all(value is False for value in self.theorem["legacy_physical_roles_not_transferred"].values()))
        self.assertFalse(self.theorem["qw2191_discharged_by_this_theorem"])
        self.assertFalse(self.theorem["toe_closure_claimed"])
        self.assertFalse(self.theorem["apd_source_exported"])

    def test_gatekeepers_and_docs(self) -> None:
        self.assertTrue(all(self.payload["gatekeeper_checks"].values()))
        self.assertEqual(self.payload["rg_non_duplication_audit"]["tool"], "rg")
        self.assertIn("P2608/S1558", MD.read_text(encoding="utf-8"))
        self.assertIn("P2608/S1558", (ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md").read_text(encoding="utf-8"))
        self.assertIn("P2608/S1558", (ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md").read_text(encoding="utf-8"))


if __name__ == "__main__":
    unittest.main()
