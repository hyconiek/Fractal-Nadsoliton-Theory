from __future__ import annotations

import json
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT))

from p2610_s1560_p2601_p2608_critical_revalidation_audit import MD, OUT, append_doc_sections, build_payload, write_markdown


class P2610P2601P2608CriticalRevalidationAuditTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.payload = build_payload()
        OUT.write_text(json.dumps(cls.payload, indent=2, sort_keys=True), encoding="utf-8")
        write_markdown(cls.payload)
        append_doc_sections()
        cls.theorem = cls.payload["p2601_p2608_critical_revalidation_audit"]["theorem_export"]

    def test_identity_and_scope(self) -> None:
        self.assertEqual(self.payload["packet_id"], "P2610")
        self.assertEqual(self.payload["stage_id"], "S1560")
        self.assertIn("CRITICAL_REVALIDATION", self.payload["status"])
        self.assertEqual(len(self.theorem["review_rows"]), 8)

    def test_critique_quarantines_weak_exports(self) -> None:
        quarantined = set(self.theorem["quarantined_packet_ids_after_revalidation"])
        self.assertTrue({"P2601", "P2602", "P2605", "P2607", "P2608"}.issubset(quarantined))
        self.assertIn("P2607_GF2_MATRIX_HAS_PHYSICAL_ORIGIN", self.theorem["strong_export_policy_after_p2610"]["legacy_to_strict_completion_bridge"])
        self.assertIn("ROLE_SEMANTICS", self.theorem["strong_export_policy_after_p2610"]["role_bearing_ltotal"])

    def test_retains_only_scoped_items(self) -> None:
        accepted = set(self.theorem["accepted_packet_ids_after_revalidation"])
        self.assertIn("P2603", accepted)
        self.assertIn("P2604", accepted)
        self.assertIn("P2606", accepted)
        self.assertEqual(self.theorem["operational_overrides"]["P2606"], "RETAIN_AS_NUMERICAL_RESIDUAL_COMPONENT_ONLY")

    def test_no_new_exports(self) -> None:
        self.assertFalse(self.theorem["strict_damping_beta_eta_source_newly_exported"])
        self.assertFalse(self.theorem["bridge_completion_newly_exported"])
        self.assertFalse(self.theorem["role_bearing_ltotal_newly_exported"])
        self.assertFalse(self.theorem["role_transfer_theorem_newly_exported"])
        self.assertFalse(self.theorem["legacy_physical_role_transfer_exported"])
        self.assertFalse(self.theorem["qw2191_discharged_by_this_audit"])
        self.assertFalse(self.theorem["toe_closure_claimed"])
        self.assertFalse(self.theorem["apd_source_exported"])

    def test_gatekeepers_and_docs(self) -> None:
        self.assertTrue(all(self.payload["gatekeeper_checks"].values()))
        self.assertEqual(self.payload["rg_non_duplication_audit"]["tool"], "rg")
        self.assertIn("P2610/S1560", MD.read_text(encoding="utf-8"))
        self.assertIn("P2610/S1560", (ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md").read_text(encoding="utf-8"))
        self.assertIn("P2610/S1560", (ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md").read_text(encoding="utf-8"))


if __name__ == "__main__":
    unittest.main()
