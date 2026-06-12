from __future__ import annotations

import json
import subprocess
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent
SCRIPT = ROOT / "p2643_s1593_inverse_hierarchy_beta_threshold_role_rejection_certificate.py"
OUT = ROOT / "generated" / "p2643_s1593_inverse_hierarchy_beta_threshold_role_rejection_certificate.json"
MD = ROOT / "generated" / "p2643_s1593_inverse_hierarchy_beta_threshold_role_rejection_certificate.md"


class P2643InverseHierarchyBetaThresholdRoleRejectionTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        subprocess.run([sys.executable, str(SCRIPT)], cwd=ROOT, check=True)
        cls.payload = json.loads(OUT.read_text(encoding="utf-8"))

    def test_content_first_grep_is_present(self) -> None:
        audit = self.payload["semantic_rg_antiduplication_audit"]
        self.assertIn("content-first", audit["mode"])
        self.assertIn("inverse_hierarchy_role_content", audit["patterns"])
        self.assertIn("beta_threshold_source_content", audit["patterns"])
        self.assertTrue(all(data["count"] > 0 for data in audit["patterns"].values()))

    def test_beta_threshold_theorem_is_computed(self) -> None:
        theorem = self.payload["beta_threshold_theorem"]
        self.assertIn("R(beta)", theorem["ratio_formula"])
        self.assertGreater(theorem["phase_ratio_C"], 3.7)
        self.assertGreater(theorem["far_power_7_eta"], 33.0)
        self.assertGreater(theorem["beta_critical_exact_role_boundary"], 0.09)
        self.assertLess(theorem["beta_critical_exact_role_boundary"], 0.10)
        self.assertEqual(theorem["derivative_sign_on_beta_nonnegative"], "negative")

    def test_ratio_table_separates_legacy_strict_and_micro_beta(self) -> None:
        rows = {row["label"]: row for row in self.payload["ratio_table"]}
        self.assertTrue(rows["legacy_beta_tors"]["preserves_unchanged_inverse_hierarchy_proxy"])
        self.assertFalse(rows["strict_beta"]["preserves_unchanged_inverse_hierarchy_proxy"])
        self.assertAlmostEqual(rows["strict_beta"]["ratio_abs_k7_over_k1"], 0.2182298586, places=9)
        self.assertFalse(rows["micro_beta_median"]["preserves_unchanged_inverse_hierarchy_proxy"])
        self.assertGreater(rows["micro_beta_median"]["beta"], rows["strict_beta"]["beta"])

    def test_role_verdict_rejects_unchanged_transfer_without_closure(self) -> None:
        verdict = self.payload["role_verdict"]
        self.assertFalse(verdict["gates"]["strict_beta_below_critical"])
        self.assertFalse(verdict["gates"]["strict_beta_preserves_ratio"])
        self.assertFalse(verdict["gates"]["micro_beta_below_critical"])
        self.assertIn("REJECT_UNCHANGED_TRANSFER", verdict["unchanged_inverse_hierarchy_role_status"])
        self.assertFalse(verdict["full_kernel_now"])
        self.assertFalse(verdict["toe_closure_now"])
        self.assertTrue(all(value is False for value in self.payload["negative_export_flags"].values()))

    def test_docs_and_recommendation_are_updated(self) -> None:
        self.assertIn("modified/compressed successor", self.payload["role_verdict"]["next_honest_step"])
        self.assertIn("P2643/S1593", MD.read_text(encoding="utf-8"))
        self.assertIn("P2643/S1593", (ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md").read_text(encoding="utf-8"))
        self.assertIn("P2643/S1593", (ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md").read_text(encoding="utf-8"))


if __name__ == "__main__":
    unittest.main()
