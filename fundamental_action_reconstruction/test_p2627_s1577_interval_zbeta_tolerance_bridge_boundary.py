from __future__ import annotations

import json
import subprocess
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent
SCRIPT = ROOT / "p2627_s1577_interval_zbeta_tolerance_bridge_boundary.py"
OUT = ROOT / "generated" / "p2627_s1577_interval_zbeta_tolerance_bridge_boundary.json"
MD = ROOT / "generated" / "p2627_s1577_interval_zbeta_tolerance_bridge_boundary.md"


class P2627IntervalZbetaToleranceBridgeBoundaryTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        subprocess.run([sys.executable, str(SCRIPT)], cwd=ROOT, check=True)
        cls.payload = json.loads(OUT.read_text(encoding="utf-8"))

    def test_content_first_grep_is_present(self) -> None:
        audit = self.payload["semantic_rg_antiduplication_audit"]
        self.assertIn("content-first", audit["mode"])
        self.assertTrue(all(data["count"] > 0 for data in audit["patterns"].values()))

    def test_tolerance_formula_has_expected_strict_interval(self) -> None:
        strict_interval = self.payload["tolerance_certificate"]["strict_epsilon_interval"]
        self.assertEqual(strict_interval["domain"], "0 < d <= 10")
        self.assertLess(strict_interval["lower"], 100.0)
        self.assertGreater(strict_interval["upper"], 100.0)
        self.assertLess(strict_interval["upper"] - strict_interval["lower"], 3.0)
        self.assertIn("maximized at D", strict_interval["proof"])

    def test_qw2064_median_only_passes_loose_practical_tolerance(self) -> None:
        rows = {row["source"]: row for row in self.payload["tolerance_certificate"]["rows"]}
        self.assertEqual(rows["target"]["relative_error_at_dmax"], 0.0)
        self.assertFalse(rows["micro_global_median"]["passes_strict_1_percent_tolerance"])
        self.assertTrue(rows["micro_global_median"]["passes_practical_15_percent_tolerance"])
        self.assertFalse(rows["bin_q25"]["passes_practical_15_percent_tolerance"])
        self.assertFalse(rows["bin_q75"]["passes_practical_15_percent_tolerance"])

    def test_reported_interval_not_admitted_as_strict_or_practical_bridge(self) -> None:
        cert = self.payload["tolerance_certificate"]
        self.assertFalse(cert["reported_interval_subset_of_strict_1_percent_admission"])
        self.assertFalse(cert["reported_interval_subset_of_practical_15_percent_admission"])

    def test_no_source_promotion_and_recommendation(self) -> None:
        acceptance = self.payload["source_acceptance"]
        self.assertFalse(acceptance["accepts_positive_beta_renormalization_source"])
        self.assertTrue(acceptance["accepts_interval_valued_support_only"])
        self.assertIn("target_independent_of_selected_kernel", acceptance["failed_gates"])
        self.assertIn("no_wide_ci_warning", acceptance["failed_gates"])
        self.assertTrue(all(value is False for value in self.payload["negative_export_flags"].values()))
        self.assertIn("approximate effective-kernel theorem", self.payload["current_blockers_and_recommendation"]["recommended_next_honest_step"])

    def test_docs_updated(self) -> None:
        self.assertIn("P2627/S1577", MD.read_text(encoding="utf-8"))
        self.assertIn("P2627/S1577", (ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md").read_text(encoding="utf-8"))
        self.assertIn("P2627/S1577", (ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md").read_text(encoding="utf-8"))


if __name__ == "__main__":
    unittest.main()
