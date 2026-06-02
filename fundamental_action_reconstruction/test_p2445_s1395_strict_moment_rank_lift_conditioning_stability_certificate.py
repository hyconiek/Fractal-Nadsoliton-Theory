import json
import subprocess
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent
SCRIPT = ROOT / "p2445_s1395_strict_moment_rank_lift_conditioning_stability_certificate.py"
OUT = ROOT / "generated" / "p2445_s1395_strict_moment_rank_lift_conditioning_stability_certificate.json"
MD = ROOT / "generated" / "p2445_s1395_strict_moment_rank_lift_conditioning_stability_certificate.md"
P2444 = ROOT / "generated" / "p2444_s1394_strict_moment_rank_lift_conditioning_certificate.json"


class P2445StrictMomentRankLiftConditioningStabilityCertificateTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        if not P2444.exists():
            subprocess.run([sys.executable, str(ROOT / "p2444_s1394_strict_moment_rank_lift_conditioning_certificate.py")], check=True)
        subprocess.run([sys.executable, str(SCRIPT)], check=True)
        cls.payload = json.loads(OUT.read_text(encoding="utf-8"))
        cls.cert = cls.payload["strict_moment_rank_lift_conditioning_stability_certificate"]
        cls.theorem = cls.cert["theorem_export"]

    def test_identity_outputs_and_grid_size(self) -> None:
        self.assertEqual(self.payload["packet_id"], "P2445")
        self.assertEqual(self.payload["stage_id"], "S1395")
        self.assertEqual(self.payload["status"], "PASS_STRICT_MOMENT_RANK_LIFT_CONDITIONING_STABILITY_NO_SOURCE_THEOREM")
        self.assertTrue(OUT.exists())
        self.assertTrue(MD.exists())
        self.assertEqual(self.theorem["configuration_count"], 12)
        self.assertEqual(self.theorem["finite_difference_step_scales"], [0.5, 1.0, 2.0, 5.0])
        self.assertEqual(self.theorem["quadrature_steps"], [10000, 20000, 40000])

    def test_conditioning_stability_frontier(self) -> None:
        self.assertEqual(self.theorem["baseline_best_candidate_id"], "K_at_d_1")
        self.assertTrue(self.theorem["all_configurations_preserve_best_candidate"])
        self.assertTrue(self.theorem["all_configurations_preserve_robust_set"])
        self.assertTrue(self.theorem["all_configurations_preserve_top_six_order"])
        self.assertEqual(self.theorem["baseline_robust_candidate_ids"], ["K_at_d_1", "K_at_d_0", "K_at_d_0.5", "K_at_d_2", "M0", "K_at_d_0.25"])
        self.assertGreater(
            self.theorem["minimum_robust_candidate_volume_across_grid"],
            self.theorem["conditioning_threshold_normalized_volume"],
        )
        self.assertLess(
            self.theorem["maximum_nonrobust_candidate_volume_across_grid"],
            self.theorem["conditioning_threshold_normalized_volume"],
        )

    def test_each_configuration_preserves_expected_rows(self) -> None:
        for result in self.theorem["configuration_results"]:
            self.assertEqual(result["best_candidate_id"], "K_at_d_1")
            self.assertEqual(result["robust_candidate_count"], 6)
            self.assertEqual(result["weakest_candidate_id"], "K_at_d_20")
            self.assertEqual(result["top_six_candidate_ids"], self.theorem["baseline_top_six_candidate_ids"])

    def test_hard_limits(self) -> None:
        self.assertFalse(self.theorem["strict_observable_source_constraint_exported_by_this_certificate"])
        self.assertFalse(self.theorem["gauge_fixing_theorem_exported_by_this_certificate"])
        self.assertFalse(self.theorem["strict_physical_value_generator_exported"])
        self.assertFalse(self.theorem["qw2191_discharged"])
        self.assertFalse(self.theorem["role_bearing_ltotal_exported"])
        self.assertFalse(self.theorem["toe_closure_exported"])
        self.assertTrue(all(self.payload["gatekeeper_checks"].values()))

    def test_docs(self) -> None:
        self.assertIn("conditioning stability certificate", MD.read_text(encoding="utf-8"))
        self.assertIn("P2445/S1395", (ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md").read_text(encoding="utf-8"))
        self.assertIn("P2445/S1395", (ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md").read_text(encoding="utf-8"))


if __name__ == "__main__":
    unittest.main()
