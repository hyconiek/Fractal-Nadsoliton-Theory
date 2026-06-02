import json
import subprocess
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent
SCRIPT = ROOT / "p2448_s1398_strict_pointwise_rank_lift_global_stationary_census_certificate.py"
OUT = ROOT / "generated" / "p2448_s1398_strict_pointwise_rank_lift_global_stationary_census_certificate.json"
MD = ROOT / "generated" / "p2448_s1398_strict_pointwise_rank_lift_global_stationary_census_certificate.md"
P2447 = ROOT / "generated" / "p2447_s1397_strict_pointwise_rank_lift_stationary_refinement_certificate.json"


class P2448StrictPointwiseRankLiftGlobalStationaryCensusCertificateTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        if not P2447.exists():
            subprocess.run([sys.executable, str(ROOT / "p2447_s1397_strict_pointwise_rank_lift_stationary_refinement_certificate.py")], check=True)
        subprocess.run([sys.executable, str(SCRIPT)], check=True)
        cls.payload = json.loads(OUT.read_text(encoding="utf-8"))
        cls.cert = cls.payload["strict_pointwise_rank_lift_global_stationary_census_certificate"]
        cls.theorem = cls.cert["theorem_export"]
        cls.census = cls.theorem["census"]

    def test_identity_outputs_and_census_size(self) -> None:
        self.assertEqual(self.payload["packet_id"], "P2448")
        self.assertEqual(self.payload["stage_id"], "S1398")
        self.assertEqual(self.payload["status"], "PASS_STRICT_POINTWISE_GLOBAL_STATIONARY_CENSUS_SELECTOR_OBSTRUCTION_NO_SOURCE_THEOREM")
        self.assertTrue(OUT.exists())
        self.assertTrue(MD.exists())
        self.assertEqual(self.census["scan_point_count"], 5001)
        self.assertEqual(self.census["sign_change_bracket_count"], 3)

    def test_stationary_census_matches_p2447(self) -> None:
        self.assertEqual(len(self.census["stationary_roots"]), 3)
        self.assertEqual(self.census["local_maximum_count"], 1)
        self.assertEqual(self.census["local_minimum_count"], 2)
        best = self.census["best_stationary_or_boundary_row"]
        self.assertEqual(best["kind"], "local_maximum")
        self.assertGreater(best["d"], 0.785)
        self.assertLess(best["d"], 0.786)
        self.assertTrue(self.theorem["best_matches_p2447_refinement"])
        self.assertTrue(self.theorem["best_volume_matches_p2447_refinement"])
        self.assertTrue(self.theorem["global_best_dominates_boundaries"])
        self.assertGreater(self.theorem["global_best_gap_over_nonbest_stationary_or_boundary_ceiling"], 0.002)
        self.assertTrue(self.theorem["all_stationary_derivatives_small"])
        self.assertTrue(self.theorem["all_stationary_roots_classified"])

    def test_hard_limits(self) -> None:
        self.assertFalse(self.theorem["analytic_interval_root_exclusion_theorem_exported_by_this_certificate"])
        self.assertFalse(self.theorem["pointwise_coordinate_selector_exported_by_this_certificate"])
        self.assertFalse(self.theorem["strict_observable_source_constraint_exported_by_this_certificate"])
        self.assertFalse(self.theorem["gauge_slice_theorem_exported_by_this_certificate"])
        self.assertFalse(self.theorem["strict_physical_value_generator_exported"])
        self.assertFalse(self.theorem["qw2191_discharged"])
        self.assertFalse(self.theorem["role_bearing_ltotal_exported"])
        self.assertFalse(self.theorem["toe_closure_exported"])
        self.assertTrue(all(self.payload["gatekeeper_checks"].values()))

    def test_docs(self) -> None:
        self.assertIn("global stationary-census certificate", MD.read_text(encoding="utf-8"))
        self.assertIn("P2448/S1398", (ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md").read_text(encoding="utf-8"))
        self.assertIn("P2448/S1398", (ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md").read_text(encoding="utf-8"))


if __name__ == "__main__":
    unittest.main()
