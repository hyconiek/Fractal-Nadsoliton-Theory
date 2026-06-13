from __future__ import annotations

import json
import subprocess
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent
SCRIPT = ROOT / "p2703_s1653_release_8_1_and_9_3s_selector_claim_reconciliation_audit.py"
OUT = ROOT / "generated" / "p2703_s1653_release_8_1_and_9_3s_selector_claim_reconciliation_audit.json"
MD = ROOT / "generated" / "p2703_s1653_release_8_1_and_9_3s_selector_claim_reconciliation_audit.md"


class P2703ReleaseSelectorClaimReconciliationAuditTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        subprocess.run([sys.executable, str(SCRIPT)], cwd=ROOT, check=True)
        cls.payload = json.loads(OUT.read_text(encoding="utf-8"))

    def test_release_8_1_positive_claims_are_detected(self) -> None:
        flags = self.payload["release_8_1_findings"]["claim_flags"]
        self.assertTrue(flags["internal_selector_source_claim"])
        self.assertTrue(flags["p1343_claim"])
        self.assertTrue(flags["p1348_global_closure_claim"])
        self.assertTrue(flags["no_false_pass_limit"])

    def test_release_9_3s_not_silently_promoted(self) -> None:
        search = self.payload["release_search"]
        self.assertIn("rg", search["tooling"])
        self.assertFalse(search["release_9_3s_document_found"])
        self.assertTrue(search["r9_theorem_draft_present"])
        r9 = self.payload["r9_related_findings"]
        self.assertTrue(r9["p1293_present"])
        self.assertTrue(r9["p1293_status_draft"])
        self.assertTrue(r9["p1293_closure_policy_blocks_closure"])

    def test_matrix_preserves_current_blocks(self) -> None:
        matrix = self.payload["reconciliation_matrix"]
        self.assertEqual(len(matrix), 5)
        self.assertTrue(any(row["older_support_real"] for row in matrix))
        self.assertTrue(all(row["current_block_removed"] is False for row in matrix))
        decision = self.payload["decision"]
        self.assertFalse(decision["does_older_release_material_remove_current_blocks"])
        self.assertIn("P2704", decision["next_honest_step"])
        self.assertTrue(all(value is False for value in self.payload["negative_export_flags"].values()))

    def test_docs_updated(self) -> None:
        self.assertIn("P2703/S1653", MD.read_text(encoding="utf-8"))
        self.assertIn("P2703/S1653", (ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md").read_text(encoding="utf-8"))
        self.assertIn("P2703/S1653", (ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md").read_text(encoding="utf-8"))
        self.assertIn("P2703/S1653", (ROOT.parent / "AGENTS.md").read_text(encoding="utf-8"))


if __name__ == "__main__":
    unittest.main()
