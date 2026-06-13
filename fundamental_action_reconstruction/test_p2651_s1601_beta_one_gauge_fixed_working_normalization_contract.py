from __future__ import annotations

import json
import subprocess
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent
SCRIPT = ROOT / "p2651_s1601_beta_one_gauge_fixed_working_normalization_contract.py"
OUT = ROOT / "generated" / "p2651_s1601_beta_one_gauge_fixed_working_normalization_contract.json"
MD = ROOT / "generated" / "p2651_s1601_beta_one_gauge_fixed_working_normalization_contract.md"


class P2651BetaOneGaugeFixedWorkingNormalizationContractTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        subprocess.run([sys.executable, str(SCRIPT)], cwd=ROOT, check=True)
        cls.payload = json.loads(OUT.read_text(encoding="utf-8"))

    def test_content_first_grep_and_upstream_consistency(self) -> None:
        audit = self.payload["semantic_rg_antiduplication_audit"]
        self.assertIn("content-first", audit["mode"])
        self.assertIn("gauge_fixed_beta_content", audit["patterns"])
        self.assertIn("rescaling_invariant_content", audit["patterns"])
        self.assertIn("role_boundary_content", audit["patterns"])
        self.assertTrue(all(data["count"] > 0 for data in audit["patterns"].values()))
        upstream = self.payload["upstream_consistency"]
        self.assertTrue(upstream["p2649_no_beta_source_routes_pass"])
        self.assertTrue(upstream["p2650_no_canonical_unit_candidates_pass"])
        self.assertTrue(upstream["p2650_recommends_gauge_fixed_working_normalization"])
        self.assertTrue(upstream["p2647_harness_ready_but_no_real_holdout"])
        self.assertTrue(upstream["p2648_margin_rule_ready_but_no_blind_data"])

    def test_orbit_invariance_and_tail_ratio_warning(self) -> None:
        orbit = self.payload["orbit_invariance_witness"]
        self.assertTrue(orbit["all_rows_invariant_to_roundoff"])
        self.assertLess(orbit["max_envelope_error"], 1e-15)
        tail = self.payload["tail_ratio_gauge_witness"]
        self.assertLess(tail["max_gauge_respecting_error"], 1e-15)
        self.assertTrue(tail["raw_substitution_has_nonzero_errors"])
        self.assertGreater(tail["max_raw_distance_substitution_error"], 0.1)

    def test_claim_boundary_contract(self) -> None:
        rows = {row["claim"]: row for row in self.payload["claim_boundary_matrix"]}
        self.assertTrue(rows["strict_beta_equals_one_working_gauge"]["allowed_under_contract"])
        self.assertTrue(rows["modified_compressed_successor_semantics"]["allowed_under_contract"])
        self.assertTrue(rows["blind_holdout_empirical_support"]["allowed_under_contract"])
        self.assertFalse(rows["target_independent_beta_source"]["allowed_under_contract"])
        self.assertFalse(rows["unchanged_legacy_inverse_hierarchy_transfer"]["allowed_under_contract"])
        self.assertFalse(rows["role_bearing_ltotal_damping_term"]["allowed_under_contract"])

    def test_no_source_or_ltotal_closure_and_docs_updated(self) -> None:
        decision = self.payload["closure_decision"]
        self.assertTrue(decision["beta_one_working_gauge_allowed"])
        self.assertFalse(decision["beta_source_exported_now"])
        self.assertFalse(decision["role_bearing_ltotal_now"])
        self.assertFalse(decision["toe_closure_now"])
        self.assertTrue(all(value is False for value in self.payload["negative_export_flags"].values()))
        self.assertIn("P2651/S1601", MD.read_text(encoding="utf-8"))
        self.assertIn("P2651/S1601", (ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md").read_text(encoding="utf-8"))
        self.assertIn("P2651/S1601", (ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md").read_text(encoding="utf-8"))


if __name__ == "__main__":
    unittest.main()
