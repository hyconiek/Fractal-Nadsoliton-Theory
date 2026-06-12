from __future__ import annotations

import json
import subprocess
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent
SCRIPT = ROOT / "p2648_s1598_frozen_kernel_holdout_statistical_margin_power_certificate.py"
OUT = ROOT / "generated" / "p2648_s1598_frozen_kernel_holdout_statistical_margin_power_certificate.json"
MD = ROOT / "generated" / "p2648_s1598_frozen_kernel_holdout_statistical_margin_power_certificate.md"


class P2648FrozenKernelHoldoutStatisticalMarginPowerTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        subprocess.run([sys.executable, str(SCRIPT)], cwd=ROOT, check=True)
        cls.payload = json.loads(OUT.read_text(encoding="utf-8"))

    def test_content_first_grep_and_upstream_are_present(self) -> None:
        audit = self.payload["semantic_rg_antiduplication_audit"]
        self.assertIn("content-first", audit["mode"])
        self.assertIn("statistical_margin_power_content", audit["patterns"])
        self.assertIn("frozen_holdout_decision_content", audit["patterns"])
        self.assertIn("legacy_strict_control_content", audit["patterns"])
        self.assertTrue(all(data["count"] > 0 for data in audit["patterns"].values()))
        upstream = self.payload["upstream_consistency"]
        self.assertTrue(upstream["p2646_has_preregistered_rows"])
        self.assertTrue(upstream["p2647_harness_ready_no_real_holdout"])
        self.assertTrue(upstream["p2647_firewall_passes"])
        self.assertTrue(upstream["p2647_requires_standard_errors"])

    def test_familywise_rule_and_margin_budget(self) -> None:
        cert = self.payload["statistical_margin_certificate"]
        rule = cert["decision_rule"]
        self.assertEqual(rule["total_inequalities"], 8)
        self.assertGreater(rule["bonferroni_one_sided_z"], 2.4)
        self.assertLess(rule["bonferroni_one_sided_z"], 2.6)
        self.assertIn("z*ratio_standard_error", rule["ratio_pass_rule_with_uncertainty"])
        self.assertIn("z*slope_standard_error", rule["slope_pass_rule_with_uncertainty"])
        self.assertTrue(cert["strict_all_nominal_margins_positive"])
        self.assertTrue(cert["legacy_all_nominal_margins_negative"])
        self.assertGreater(cert["minimum_max_ratio_standard_error_for_strict_pass"], 0.0)
        self.assertGreater(cert["minimum_max_slope_standard_error_for_strict_pass"], 0.0)

    def test_pairwise_rows_keep_legacy_on_wrong_side(self) -> None:
        rows = self.payload["statistical_margin_certificate"]["rows"]
        self.assertEqual(len(rows), 4)
        for row in rows:
            self.assertGreater(row["strict_ratio_margin"], 0.0)
            self.assertGreater(row["strict_slope_margin"], 0.0)
            self.assertLess(row["legacy_ratio_margin"], 0.0)
            self.assertLess(row["legacy_slope_margin"], 0.0)
            self.assertTrue(row["legacy_ratio_is_already_wrong_side"])
            self.assertTrue(row["legacy_slope_is_already_wrong_side"])
            self.assertTrue(row["midpoint_fails_strict_inequality_before_uncertainty"])

    def test_no_empirical_or_source_closure_and_docs_updated(self) -> None:
        decision = self.payload["closure_decision"]
        self.assertIn("STATISTICAL_MARGIN_RULE_READY", decision["decision"])
        self.assertFalse(decision["real_blind_holdout_executed"])
        self.assertFalse(decision["full_kernel_now"])
        self.assertFalse(decision["toe_closure_now"])
        self.assertTrue(decision["can_upgrade_p2647_harness"])
        self.assertTrue(all(value is False for value in self.payload["negative_export_flags"].values()))
        self.assertIn("P2648/S1598", MD.read_text(encoding="utf-8"))
        self.assertIn("P2648/S1598", (ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md").read_text(encoding="utf-8"))
        self.assertIn("P2648/S1598", (ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md").read_text(encoding="utf-8"))


if __name__ == "__main__":
    unittest.main()
