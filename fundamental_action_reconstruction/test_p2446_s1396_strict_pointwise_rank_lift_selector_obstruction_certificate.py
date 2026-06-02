import json
import subprocess
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent
SCRIPT = ROOT / "p2446_s1396_strict_pointwise_rank_lift_selector_obstruction_certificate.py"
OUT = ROOT / "generated" / "p2446_s1396_strict_pointwise_rank_lift_selector_obstruction_certificate.json"
MD = ROOT / "generated" / "p2446_s1396_strict_pointwise_rank_lift_selector_obstruction_certificate.md"
P2445 = ROOT / "generated" / "p2445_s1395_strict_moment_rank_lift_conditioning_stability_certificate.json"


class P2446StrictPointwiseRankLiftSelectorObstructionCertificateTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        if not P2445.exists():
            subprocess.run([sys.executable, str(ROOT / "p2445_s1395_strict_moment_rank_lift_conditioning_stability_certificate.py")], check=True)
        subprocess.run([sys.executable, str(SCRIPT)], check=True)
        cls.payload = json.loads(OUT.read_text(encoding="utf-8"))
        cls.cert = cls.payload["strict_pointwise_rank_lift_selector_obstruction_certificate"]
        cls.theorem = cls.cert["theorem_export"]

    def test_identity_outputs_and_scan_size(self) -> None:
        self.assertEqual(self.payload["packet_id"], "P2446")
        self.assertEqual(self.payload["stage_id"], "S1396")
        self.assertEqual(self.payload["status"], "PASS_STRICT_POINTWISE_RANK_LIFT_SELECTOR_OBSTRUCTION_NO_SOURCE_THEOREM")
        self.assertTrue(OUT.exists())
        self.assertTrue(MD.exists())
        self.assertEqual(self.theorem["scan_point_count"], 2001)
        self.assertEqual(self.theorem["inherited_p2445_best_candidate_id"], "K_at_d_1")

    def test_pointwise_selector_obstruction(self) -> None:
        self.assertTrue(self.theorem["d1_is_robust"])
        self.assertFalse(self.theorem["d1_is_global_conditioning_maximum_on_scan"])
        self.assertFalse(self.theorem["d1_is_local_conditioning_maximum_on_window"])
        self.assertGreater(self.theorem["best_volume_minus_d1_volume"], 0.0)
        self.assertGreater(self.theorem["local_best_volume_minus_d1_volume"], 0.0)
        self.assertAlmostEqual(self.theorem["global_best_pointwise_row"]["d"], 0.785, places=12)
        self.assertEqual(self.theorem["local_best_pointwise_row"]["d"], self.theorem["global_best_pointwise_row"]["d"])
        self.assertGreater(self.theorem["robust_pointwise_count_in_local_window"], 1)

    def test_hard_limits(self) -> None:
        self.assertFalse(self.theorem["pointwise_coordinate_selector_exported_by_this_certificate"])
        self.assertFalse(self.theorem["strict_observable_source_constraint_exported_by_this_certificate"])
        self.assertFalse(self.theorem["gauge_fixing_theorem_exported_by_this_certificate"])
        self.assertFalse(self.theorem["strict_physical_value_generator_exported"])
        self.assertFalse(self.theorem["qw2191_discharged"])
        self.assertFalse(self.theorem["role_bearing_ltotal_exported"])
        self.assertFalse(self.theorem["toe_closure_exported"])
        self.assertTrue(all(self.payload["gatekeeper_checks"].values()))

    def test_docs(self) -> None:
        self.assertIn("selector obstruction certificate", MD.read_text(encoding="utf-8"))
        self.assertIn("P2446/S1396", (ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md").read_text(encoding="utf-8"))
        self.assertIn("P2446/S1396", (ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md").read_text(encoding="utf-8"))


if __name__ == "__main__":
    unittest.main()
