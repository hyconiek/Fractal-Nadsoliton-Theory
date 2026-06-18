import json
import subprocess
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
P2865_SCRIPT = ROOT / "p2865_s1815_singleton_localizer_coefficient_product_obligation_no_go_audit.py"
SCRIPT = ROOT / "p2866_s1816_product_theorem_candidate_pair_no_source_audit.py"
JSON_PATH = ROOT / "generated" / "p2866_s1816_product_theorem_candidate_pair_no_source_audit.json"
MD_PATH = ROOT / "generated" / "p2866_s1816_product_theorem_candidate_pair_no_source_audit.md"


class P2866ProductTheoremCandidatePairNoSourceAuditTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        for script in (P2865_SCRIPT, SCRIPT):
            subprocess.run([sys.executable, str(script)], check=True, cwd=ROOT, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
        cls.payload = json.loads(JSON_PATH.read_text(encoding="utf-8"))
        cls.audit = cls.payload["product_theorem_candidate_pair_no_source_audit"]

    def test_status_and_p2865_input(self):
        self.assertEqual(
            self.payload["status"],
            "P2866_PRODUCT_THEOREM_CANDIDATE_PAIR_NO_SOURCE_AUDIT_NO_CLOSURE",
        )
        self.assertEqual(
            self.audit["input_status_rechecked"],
            "P2865_SINGLETON_LOCALIZER_COEFFICIENT_PRODUCT_OBLIGATION_NO_GO_AUDIT_NO_CLOSURE",
        )

    def test_nine_candidate_pairs_and_unique_exact_imported_pair(self):
        rows = self.audit["candidate_pair_matrix"]
        self.assertEqual(len(rows), 9)
        exact = self.audit["exact_representation_rows"]
        self.assertEqual(len(exact), 1)
        self.assertEqual(exact[0]["localizer"], "singleton_d11_localizer")
        self.assertEqual(exact[0]["coefficient_law"], "imported_prime5_coefficient_9_over_5")
        self.assertIn("localizer_not_strict_sourced", exact[0]["failure_modes"])
        self.assertIn("coefficient_not_strict_sourced", exact[0]["failure_modes"])

    def test_sourced_pairs_fail_exactness(self):
        sourced_side_rows = [
            row
            for row in self.audit["candidate_pair_matrix"]
            if row["localizer_strict_sourced_nonpremise"] or row["coefficient_strict_sourced_nonpremise"]
        ]
        self.assertTrue(sourced_side_rows)
        exact_and_fully_sourced = [
            row
            for row in sourced_side_rows
            if row["exact_representation_of_target"]
            and row["localizer_strict_sourced_nonpremise"]
            and row["coefficient_strict_sourced_nonpremise"]
        ]
        self.assertEqual(exact_and_fully_sourced, [])

    def test_no_false_exports(self):
        self.assertEqual(self.audit["accepted_candidate_count"], 0)
        self.assertFalse(self.payload["acceptance_matrix"]["exports_boundary_source_law"])
        self.assertFalse(self.payload["acceptance_matrix"]["exports_unit_bearing_coupling_localization_theorem"])
        flags = self.payload["decision"]["negative_export_flags"]
        self.assertFalse(flags["eta_source_exported"])
        self.assertFalse(flags["selector_closure_exported"])
        self.assertFalse(flags["ltotal_exported"])
        self.assertFalse(flags["toe_closure_exported"])

    def test_documents_updated(self):
        self.assertIn("P2866/S1816", MD_PATH.read_text(encoding="utf-8"))
        self.assertIn("P2866/S1816", (ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md").read_text(encoding="utf-8"))
        self.assertIn("P2866/S1816", (ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md").read_text(encoding="utf-8"))
        self.assertIn("P2866", (REPO / "AGENTS.md").read_text(encoding="utf-8"))


if __name__ == "__main__":
    unittest.main()
