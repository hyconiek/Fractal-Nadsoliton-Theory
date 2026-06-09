from __future__ import annotations

import json
import subprocess
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent
SCRIPT = ROOT / "p2634_s1584_strict_stability_evidence_vs_role_completeness_audit.py"
OUT = ROOT / "generated" / "p2634_s1584_strict_stability_evidence_vs_role_completeness_audit.json"
MD = ROOT / "generated" / "p2634_s1584_strict_stability_evidence_vs_role_completeness_audit.md"


class P2634StrictStabilityEvidenceVsRoleCompletenessAuditTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        subprocess.run([sys.executable, str(SCRIPT)], cwd=ROOT, check=True)
        cls.payload = json.loads(OUT.read_text(encoding="utf-8"))

    def test_content_first_grep_is_present_and_nonempty(self) -> None:
        audit = self.payload["semantic_rg_antiduplication_audit"]
        self.assertIn("content-first", audit["mode"])
        self.assertTrue(all(data["count"] > 0 for data in audit["patterns"].values()))
        self.assertIn("strict_stability_robustness_content", audit["patterns"])
        self.assertIn("legacy_role_gap_content", audit["patterns"])

    def test_stability_evidence_is_counted_as_positive_but_not_finality(self) -> None:
        ledger = self.payload["stability_evidence_ledger"]
        self.assertEqual(len(ledger), 6)
        self.assertTrue(any(row["source"] == "QW2049" and row["supports_strict_internal_stability"] for row in ledger))
        self.assertTrue(any(row["source"] == "QW2051" and row["selected_kernel_stable"] for row in ledger))
        scores = self.payload["aggregate_scores"]
        self.assertGreaterEqual(scores["strict_internal_stability_positive_classes"], 5)
        self.assertIn("QW1968", scores["fragility_warning_sources"])

    def test_role_completeness_blocks_remain_open(self) -> None:
        roles = self.payload["role_completeness_obstruction_ledger"]
        self.assertEqual(len(roles), 5)
        self.assertTrue(all(row["closed_by_stability_tests"] is False for row in roles))
        self.assertTrue(any("inverse hierarchy" in row["missing_or_guarded_role"] and row["blocking_signal"] for row in roles))
        self.assertGreaterEqual(self.payload["aggregate_scores"]["role_completeness_open_blocks"], 4)

    def test_stability_and_role_completion_are_logically_orthogonal(self) -> None:
        table = self.payload["orthogonality_truth_table"]
        self.assertEqual(len(table), 4)
        stable_incomplete = next(row for row in table if row["strict_internal_stability_evidence_passes"] and not row["legacy_to_strict_characteristic_roles_complete"])
        self.assertTrue(stable_incomplete["may_claim_working_kernel_stable"])
        self.assertFalse(stable_incomplete["may_claim_final_toe_kernel"])
        self.assertEqual(stable_incomplete["verdict"], "stable_working_successor_not_final")

    def test_acceptance_and_negative_exports(self) -> None:
        acceptance = self.payload["source_acceptance"]
        self.assertFalse(acceptance["accepts_stability_as_final_toe_kernel_completion"])
        self.assertIn("phase_frequency_node_gauge_certificate_closed", acceptance["failed_gates"])
        self.assertTrue(all(value is False for value in self.payload["negative_export_flags"].values()))

    def test_recommendation_and_docs_updated(self) -> None:
        self.assertIn("typed stability-to-role interface theorem", self.payload["recommended_next_honest_step"])
        self.assertIn("P2634/S1584", MD.read_text(encoding="utf-8"))
        self.assertIn("P2634/S1584", (ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md").read_text(encoding="utf-8"))
        self.assertIn("P2634/S1584", (ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md").read_text(encoding="utf-8"))


if __name__ == "__main__":
    unittest.main()
