from __future__ import annotations

import json
import subprocess
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent
SCRIPT = ROOT / "p2629_s1579_zbeta_normalization_gauge_obstruction.py"
OUT = ROOT / "generated" / "p2629_s1579_zbeta_normalization_gauge_obstruction.json"
MD = ROOT / "generated" / "p2629_s1579_zbeta_normalization_gauge_obstruction.md"


class P2629ZbetaNormalizationGaugeObstructionTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        subprocess.run([sys.executable, str(SCRIPT)], cwd=ROOT, check=True)
        cls.payload = json.loads(OUT.read_text(encoding="utf-8"))

    def test_content_first_grep_is_present(self) -> None:
        audit = self.payload["semantic_rg_antiduplication_audit"]
        self.assertIn("content-first", audit["mode"])
        self.assertTrue(all(data["count"] > 0 for data in audit["patterns"].values()))

    def test_invariant_ratio_is_micro_over_strict_beta(self) -> None:
        metrics = self.payload["qw_metrics"]
        self.assertAlmostEqual(metrics["micro_beta_median"], 1.1473957999384183)
        self.assertAlmostEqual(metrics["strict_beta"], 1.0)
        self.assertAlmostEqual(metrics["z_median_over_z_target"], metrics["micro_over_strict_beta"])
        self.assertGreater(metrics["relative_beta_error"], 0.14)

    def test_normalization_orbit_preserves_mismatch_ratio(self) -> None:
        orbit = self.payload["normalization_orbit_certificate"]
        invariant = orbit["invariant_ratio"]
        self.assertGreater(invariant, 1.14)
        self.assertLess(invariant, 1.15)
        for row in orbit["rows"]:
            self.assertAlmostEqual(row["median_over_target_invariant"], invariant)
        self.assertIn("beta_uv -> lambda beta_uv", orbit["proof"])

    def test_no_exact_source_promotion(self) -> None:
        gate = self.payload["exact_source_gate"]
        self.assertFalse(gate["accepts_uv_normalization_source"])
        self.assertFalse(gate["accepts_positive_beta_renormalization_source"])
        self.assertIn("uv_normalization_fixed_by_independent_theorem", gate["failed_gates"])
        self.assertIn("micro_beta_within_strict_1pct", gate["failed_gates"])
        self.assertTrue(all(value is False for value in self.payload["negative_export_flags"].values()))

    def test_recommendation_and_docs_updated(self) -> None:
        self.assertIn("target-independent micro normalization identity", self.payload["recommended_next_honest_step"])
        self.assertIn("P2629/S1579", MD.read_text(encoding="utf-8"))
        self.assertIn("P2629/S1579", (ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md").read_text(encoding="utf-8"))
        self.assertIn("P2629/S1579", (ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md").read_text(encoding="utf-8"))


if __name__ == "__main__":
    unittest.main()
