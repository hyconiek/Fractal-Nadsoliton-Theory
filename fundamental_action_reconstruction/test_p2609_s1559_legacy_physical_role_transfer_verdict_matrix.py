from __future__ import annotations

import json
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT))

from p2609_s1559_legacy_physical_role_transfer_verdict_matrix import MD, OUT, ROLE_CLAIMS, append_doc_sections, build_payload, write_markdown


class P2609LegacyPhysicalRoleTransferVerdictMatrixTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.payload = build_payload()
        OUT.write_text(json.dumps(cls.payload, indent=2, sort_keys=True), encoding="utf-8")
        write_markdown(cls.payload)
        append_doc_sections()
        cls.theorem = cls.payload["legacy_physical_role_transfer_verdict_matrix"]["theorem_export"]

    def test_identity_and_inheritance(self) -> None:
        self.assertEqual(self.payload["packet_id"], "P2609")
        self.assertEqual(self.payload["stage_id"], "S1559")
        self.assertIn("LEGACY_PHYSICAL_ROLE_TRANSFER_AUDIT", self.payload["status"])
        self.assertTrue(self.theorem["legacy_to_strict_completion_bridge_inherited"])
        self.assertTrue(self.theorem["strict_damping_role_transfer_inherited"])
        self.assertTrue(self.theorem["role_bearing_ltotal_for_strict_damping_inherited"])

    def test_all_legacy_claims_classified_and_rejected_now(self) -> None:
        self.assertEqual(self.theorem["legacy_claim_count"], len(ROLE_CLAIMS))
        self.assertEqual(self.theorem["accepted_legacy_claim_count"], 0)
        self.assertEqual(self.theorem["rejected_legacy_claim_count"], len(ROLE_CLAIMS))
        self.assertEqual(set(self.theorem["rejected_legacy_claim_ids"]), {claim["claim_id"] for claim in ROLE_CLAIMS})
        for row in self.theorem["legacy_claim_verdict_rows"]:
            self.assertEqual(row["verdict"], "REJECT_STRICT_TRANSFER_NOW")
            self.assertFalse(row["transfer_accepts"])
            self.assertGreater(row["missing_condition_count"], 0)
            self.assertTrue(row["required_strict_ingredient"])

    def test_scope_guards(self) -> None:
        self.assertEqual(self.theorem["survives_without_modification"], [])
        self.assertEqual(self.theorem["modified_survivals"][0]["role"], "strict_damping_compression_beta_eta_term")
        self.assertFalse(self.theorem["legacy_physical_role_transfer_exported"])
        self.assertFalse(self.theorem["legacy_electroweak_role_exported"])
        self.assertFalse(self.theorem["legacy_alpha_em_role_exported"])
        self.assertFalse(self.theorem["legacy_gravity_hierarchy_role_exported"])
        self.assertFalse(self.theorem["legacy_beta_tors_orientation_role_exported"])
        self.assertFalse(self.theorem["qw2191_discharged_by_this_packet"])
        self.assertFalse(self.theorem["toe_closure_claimed"])
        self.assertFalse(self.theorem["apd_source_exported"])

    def test_gatekeepers_and_docs(self) -> None:
        self.assertTrue(all(self.payload["gatekeeper_checks"].values()))
        self.assertEqual(self.payload["rg_non_duplication_audit"]["tool"], "rg")
        self.assertIn("P2609/S1559", MD.read_text(encoding="utf-8"))
        self.assertIn("P2609/S1559", (ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md").read_text(encoding="utf-8"))
        self.assertIn("P2609/S1559", (ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md").read_text(encoding="utf-8"))


if __name__ == "__main__":
    unittest.main()
