from __future__ import annotations

import json
import subprocess
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent
SCRIPT = ROOT / "p2650_s1600_canonical_length_uv_unit_source_candidate_exhaustion_no_go.py"
OUT = ROOT / "generated" / "p2650_s1600_canonical_length_uv_unit_source_candidate_exhaustion_no_go.json"
MD = ROOT / "generated" / "p2650_s1600_canonical_length_uv_unit_source_candidate_exhaustion_no_go.md"


class P2650CanonicalLengthUvUnitSourceCandidateExhaustionTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        subprocess.run([sys.executable, str(SCRIPT)], cwd=ROOT, check=True)
        cls.payload = json.loads(OUT.read_text(encoding="utf-8"))

    def test_content_first_grep_and_upstream_consistency(self) -> None:
        audit = self.payload["semantic_rg_antiduplication_audit"]
        self.assertIn("content-first", audit["mode"])
        self.assertIn("canonical_length_uv_content", audit["patterns"])
        self.assertIn("nadsoliton_source_content", audit["patterns"])
        self.assertIn("legacy_strict_unit_bridge_content", audit["patterns"])
        self.assertTrue(all(data["count"] > 0 for data in audit["patterns"].values()))
        upstream = self.payload["upstream_consistency"]
        self.assertTrue(upstream["p2629_normalization_gauge_unfixed"])
        self.assertGreater(upstream["p2629_micro_over_strict_ratio"], 1.0)
        self.assertGreater(upstream["p2643_beta_crit_available"], 0.0)
        self.assertTrue(upstream["p2649_no_beta_source_routes_pass"])
        self.assertTrue(upstream["p2649_demands_canonical_length_or_conservation_identity"])

    def test_all_named_canonical_unit_candidates_fail(self) -> None:
        candidates = self.payload["canonical_unit_candidates"]
        self.assertEqual(len(candidates), 7)
        self.assertFalse(any(row["passes_as_canonical_length_uv_source"] for row in candidates))
        by_candidate = {row["candidate"]: row for row in candidates}
        self.assertIn("uv_normalization_gauge_fixed", by_candidate["dimensionless_domain_unit_d_equals_1"]["missing_required_atoms"])
        self.assertIn("role_transfer_theorem", by_candidate["legacy_beta_tors_unit"]["missing_required_atoms"])
        self.assertIn("alpha_geo_to_damping_operator_source", by_candidate["alpha_geo_unit"]["missing_required_atoms"])
        self.assertIn("micro_strict_mismatch_removed", by_candidate["micro_beta_median_unit"]["missing_required_atoms"])
        self.assertIn("empirical_holdout_can_replace_source_theorem", by_candidate["empirical_tail_ratio_unit"]["missing_required_atoms"])

    def test_finite_expression_scan_remains_numerology_not_source(self) -> None:
        scan = self.payload["finite_expression_scan"]
        self.assertFalse(scan["exports_canonical_unit"])
        self.assertIn("two-factor products", scan["grammar"])
        self.assertIn("untyped", scan["verdict"])
        self.assertGreaterEqual(scan["near_rows_count"], 0)
        if scan["near_rows_count"]:
            self.assertLess(scan["minimum_distance_to_beta_one"], 0.05)

    def test_no_closure_and_docs_updated(self) -> None:
        decision = self.payload["closure_decision"]
        self.assertEqual(decision["passing_candidates"], [])
        self.assertFalse(decision["canonical_length_uv_unit_exported_now"])
        self.assertFalse(decision["beta_source_exported_now"])
        self.assertFalse(decision["full_kernel_now"])
        self.assertFalse(decision["toe_closure_now"])
        self.assertTrue(all(value is False for value in self.payload["negative_export_flags"].values()))
        self.assertIn("P2650/S1600", MD.read_text(encoding="utf-8"))
        self.assertIn("P2650/S1600", (ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md").read_text(encoding="utf-8"))
        self.assertIn("P2650/S1600", (ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md").read_text(encoding="utf-8"))


if __name__ == "__main__":
    unittest.main()
