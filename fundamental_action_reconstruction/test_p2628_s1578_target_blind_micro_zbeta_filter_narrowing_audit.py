from __future__ import annotations

import json
import subprocess
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent
SCRIPT = ROOT / "p2628_s1578_target_blind_micro_zbeta_filter_narrowing_audit.py"
OUT = ROOT / "generated" / "p2628_s1578_target_blind_micro_zbeta_filter_narrowing_audit.json"
MD = ROOT / "generated" / "p2628_s1578_target_blind_micro_zbeta_filter_narrowing_audit.md"


class P2628TargetBlindMicroZbetaFilterNarrowingAuditTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        subprocess.run([sys.executable, str(SCRIPT)], cwd=ROOT, check=True)
        cls.payload = json.loads(OUT.read_text(encoding="utf-8"))

    def test_content_first_grep_is_present(self) -> None:
        audit = self.payload["semantic_rg_antiduplication_audit"]
        self.assertIn("content-first", audit["mode"])
        self.assertTrue(all(data["count"] > 0 for data in audit["patterns"].values()))

    def test_raw_quality_summary_reads_qw2048_bins(self) -> None:
        summary = self.payload["raw_quality_summary"]
        self.assertEqual(summary["n_bins"], 17)
        self.assertLess(summary["all_bin_z_beta_min"], 20.0)
        self.assertGreater(summary["all_bin_z_beta_max"], 2000.0)
        self.assertEqual(summary["support_n_range"], [8, 85])

    def test_filter_exhaustion_has_zero_strict_acceptance(self) -> None:
        cert = self.payload["target_blind_filter_certificate"]
        self.assertGreater(cert["evaluated_filter_count"], 100)
        self.assertEqual(cert["strict_accepting_filter_count"], 0)
        self.assertIn("no filter", cert["finite_exhaustion_theorem"])

    def test_best_row_is_not_strict_even_when_median_nearish(self) -> None:
        best = self.payload["target_blind_filter_certificate"]["best_by_median_relative_error"][0]
        self.assertAlmostEqual(best["z_beta_median_of_bins"], 109.76147103809966)
        self.assertGreater(best["median_relative_error_to_100"], 0.09)
        self.assertFalse(best["accepted_as_strict_positive_beta_source"])
        self.assertFalse(best["accepted_as_practical_interval_support"])

    def test_no_source_promotion_and_next_recommendation(self) -> None:
        acceptance = self.payload["acceptance_and_recommendation"]
        self.assertFalse(acceptance["accepts_positive_beta_renormalization_source"])
        self.assertIn("micro operator", acceptance["recommended_next_honest_step"])
        self.assertTrue(all(value is False for value in self.payload["negative_export_flags"].values()))

    def test_docs_updated(self) -> None:
        self.assertIn("P2628/S1578", MD.read_text(encoding="utf-8"))
        self.assertIn("P2628/S1578", (ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md").read_text(encoding="utf-8"))
        self.assertIn("P2628/S1578", (ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md").read_text(encoding="utf-8"))


if __name__ == "__main__":
    unittest.main()
