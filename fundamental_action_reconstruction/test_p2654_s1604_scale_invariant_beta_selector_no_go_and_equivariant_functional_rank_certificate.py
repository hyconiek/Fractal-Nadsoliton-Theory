from __future__ import annotations

import json
import subprocess
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent
SCRIPT = ROOT / "p2654_s1604_scale_invariant_beta_selector_no_go_and_equivariant_functional_rank_certificate.py"
OUT = ROOT / "generated" / "p2654_s1604_scale_invariant_beta_selector_no_go_and_equivariant_functional_rank_certificate.json"
MD = ROOT / "generated" / "p2654_s1604_scale_invariant_beta_selector_no_go_and_equivariant_functional_rank_certificate.md"


class P2654ScaleInvariantBetaSelectorNoGoTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        subprocess.run([sys.executable, str(SCRIPT)], cwd=ROOT, check=True)
        cls.payload = json.loads(OUT.read_text(encoding="utf-8"))

    def test_content_first_grep_and_upstream_consistency(self) -> None:
        audit = self.payload["semantic_rg_antiduplication_audit"]
        self.assertIn("content-first", audit["mode"])
        self.assertIn("scale_invariant_selector_content", audit["patterns"])
        self.assertIn("typed_unit_source_content", audit["patterns"])
        self.assertIn("raw_anchor_warning_content", audit["patterns"])
        self.assertTrue(all(data["count"] > 0 for data in audit["patterns"].values()))
        upstream = self.payload["upstream_consistency"]
        self.assertTrue(upstream["p2649_beta_routes_blocked"])
        self.assertTrue(upstream["p2651_beta_one_only_working_gauge"])
        self.assertTrue(upstream["p2652_unit_map_source_not_exported"])
        self.assertTrue(upstream["p2653_typed_metric_uv_source_not_exported"])
        self.assertTrue(upstream["p2653_scale_orbit_equivalence_verified"])

    def test_feature_rank_no_go(self) -> None:
        no_go = self.payload["feature_rank_no_go"]
        self.assertEqual(no_go["invariant_feature_rank_on_audited_orbit"], 1)
        self.assertLess(no_go["max_invariant_difference"], 1e-12)
        self.assertTrue(no_go["raw_anchor_features_distinguish_representatives"])
        self.assertFalse(no_go["scale_invariant_selector_exists_on_audited_features"])

    def test_selector_matrix_blocks_all_current_candidates(self) -> None:
        matrix = {row["candidate"]: row for row in self.payload["selector_candidate_matrix"]}
        self.assertFalse(matrix["orbit_invariant_envelope_grid_selector"]["passes_as_target_independent_beta_source_now"])
        self.assertFalse(matrix["orbit_invariant_tail_ratio_log_slope_selector"]["passes_as_target_independent_beta_source_now"])
        self.assertTrue(matrix["raw_coordinate_envelope_anchor_selector"]["uses_external_unit_anchor"])
        self.assertFalse(matrix["raw_coordinate_envelope_anchor_selector"]["passes_as_target_independent_beta_source_now"])
        self.assertFalse(matrix["typed_metric_uv_plus_operator_identity_selector"]["passes_as_target_independent_beta_source_now"])

    def test_no_closure_and_docs_updated(self) -> None:
        decision = self.payload["closure_decision"]
        self.assertEqual(decision["passing_beta_source_candidates"], [])
        self.assertFalse(decision["scale_invariant_selector_exists_now"])
        self.assertFalse(decision["raw_anchor_promoted_to_source_now"])
        self.assertFalse(decision["canonical_unit_exported_now"])
        self.assertFalse(decision["beta_source_exported_now"])
        self.assertFalse(decision["role_bearing_ltotal_now"])
        self.assertFalse(decision["toe_closure_now"])
        self.assertTrue(all(value is False for value in self.payload["negative_export_flags"].values()))
        self.assertIn("P2654/S1604", MD.read_text(encoding="utf-8"))
        self.assertIn("P2654/S1604", (ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md").read_text(encoding="utf-8"))
        self.assertIn("P2654/S1604", (ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md").read_text(encoding="utf-8"))


if __name__ == "__main__":
    unittest.main()
