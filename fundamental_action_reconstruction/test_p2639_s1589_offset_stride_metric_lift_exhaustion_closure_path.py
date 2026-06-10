from __future__ import annotations

import json
import subprocess
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent
SCRIPT = ROOT / "p2639_s1589_offset_stride_metric_lift_exhaustion_closure_path.py"
OUT = ROOT / "generated" / "p2639_s1589_offset_stride_metric_lift_exhaustion_closure_path.json"
MD = ROOT / "generated" / "p2639_s1589_offset_stride_metric_lift_exhaustion_closure_path.md"


class P2639OffsetStrideMetricLiftExhaustionClosurePathTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        subprocess.run([sys.executable, str(SCRIPT)], cwd=ROOT, check=True)
        cls.payload = json.loads(OUT.read_text(encoding="utf-8"))

    def test_content_first_grep_is_present(self) -> None:
        audit = self.payload["semantic_rg_antiduplication_audit"]
        self.assertIn("content-first", audit["mode"])
        self.assertTrue(all(data["count"] > 0 for data in audit["patterns"].values()))
        self.assertIn("offset_stride_zero_lattice_content", audit["patterns"])
        self.assertIn("metric_coordinate_source_content", audit["patterns"])

    def test_exhaustion_finds_uv_safe_exact_lifts(self) -> None:
        exh = self.payload["offset_stride_exhaustion"]
        self.assertEqual(exh["total_candidates"], 90)
        self.assertGreater(exh["uv_safe_exact_lift_count"], 0)
        nearest = exh["nearest_local_uv_safe_exact_lift"]
        self.assertEqual(nearest["map"], "r(d)=(4/3)*d+(8/3)")
        self.assertTrue(nearest["exact_node_lift"])
        self.assertTrue(nearest["uv_domain_positive_on_grid_0_1_step_1_12"])

    def test_some_lifts_preserve_ratio_but_remain_unsourced(self) -> None:
        exh = self.payload["offset_stride_exhaustion"]
        self.assertGreater(exh["uv_safe_and_inverse_hierarchy_count"], 0)
        best = exh["best_inverse_hierarchy_uv_safe_lift"]
        self.assertTrue(best["inverse_hierarchy_ratio_above_one"])
        self.assertGreater(best["strict_lifted_abs_k7_over_k1"], 1.0)
        self.assertGreaterEqual(best["stride"], 1)

    def test_closure_decision_does_not_promote(self) -> None:
        decision = self.payload["closure_decision"]
        self.assertTrue(decision["gates"]["some_domain_safe_exact_lift_exists"])
        self.assertTrue(decision["gates"]["some_domain_safe_inverse_hierarchy_lift_exists"])
        self.assertFalse(decision["gates"]["canonical_offset_stride_source_theorem_exists"])
        self.assertFalse(decision["promote_any_lift_to_bridge_completion"])
        self.assertFalse(decision["full_kernel_now"])
        self.assertTrue(all(value is False for value in self.payload["negative_export_flags"].values()))

    def test_professorial_path_and_docs(self) -> None:
        path = self.payload["professorial_closure_path"]
        self.assertEqual(path[0]["task"], "canonical offset/stride source theorem")
        self.assertIn("topology/selector", self.payload["closure_decision"]["next_honest_step"])
        self.assertIn("P2639/S1589", MD.read_text(encoding="utf-8"))
        self.assertIn("P2639/S1589", (ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md").read_text(encoding="utf-8"))
        self.assertIn("P2639/S1589", (ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md").read_text(encoding="utf-8"))


if __name__ == "__main__":
    unittest.main()
