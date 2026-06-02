import json
import subprocess
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent
SCRIPT = ROOT / "p2447_s1397_strict_pointwise_rank_lift_stationary_refinement_certificate.py"
OUT = ROOT / "generated" / "p2447_s1397_strict_pointwise_rank_lift_stationary_refinement_certificate.json"
MD = ROOT / "generated" / "p2447_s1397_strict_pointwise_rank_lift_stationary_refinement_certificate.md"
P2446 = ROOT / "generated" / "p2446_s1396_strict_pointwise_rank_lift_selector_obstruction_certificate.json"


class P2447StrictPointwiseRankLiftStationaryRefinementCertificateTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        if not P2446.exists():
            subprocess.run([sys.executable, str(ROOT / "p2446_s1396_strict_pointwise_rank_lift_selector_obstruction_certificate.py")], check=True)
        subprocess.run([sys.executable, str(SCRIPT)], check=True)
        cls.payload = json.loads(OUT.read_text(encoding="utf-8"))
        cls.cert = cls.payload["strict_pointwise_rank_lift_stationary_refinement_certificate"]
        cls.theorem = cls.cert["theorem_export"]

    def test_identity_outputs_and_inheritance(self) -> None:
        self.assertEqual(self.payload["packet_id"], "P2447")
        self.assertEqual(self.payload["stage_id"], "S1397")
        self.assertEqual(self.payload["status"], "PASS_STRICT_POINTWISE_STATIONARY_REFINEMENT_SELECTOR_OBSTRUCTION_NO_SOURCE_THEOREM")
        self.assertTrue(OUT.exists())
        self.assertTrue(MD.exists())
        self.assertAlmostEqual(self.theorem["inherited_p2446_grid_best_d"], 0.785, places=12)

    def test_stationary_refinement(self) -> None:
        maximum = self.theorem["golden_section_maximum"]
        self.assertGreater(maximum["d"], 0.785)
        self.assertLess(maximum["d"], 0.786)
        self.assertGreater(self.theorem["conditioning_gap_over_d1"], 0.001)
        self.assertTrue(self.theorem["continuous_maximum_exceeds_d1"])
        self.assertTrue(self.theorem["continuous_maximum_inside_local_window"])
        self.assertTrue(self.theorem["continuous_maximum_exceeds_local_boundaries"])
        self.assertTrue(self.theorem["all_first_derivative_witnesses_small"])
        self.assertTrue(self.theorem["all_second_derivative_witnesses_negative"])

    def test_hard_limits(self) -> None:
        self.assertFalse(self.theorem["pointwise_coordinate_selector_exported_by_this_certificate"])
        self.assertFalse(self.theorem["strict_observable_source_constraint_exported_by_this_certificate"])
        self.assertFalse(self.theorem["gauge_slice_theorem_exported_by_this_certificate"])
        self.assertFalse(self.theorem["strict_physical_value_generator_exported"])
        self.assertFalse(self.theorem["qw2191_discharged"])
        self.assertFalse(self.theorem["role_bearing_ltotal_exported"])
        self.assertFalse(self.theorem["toe_closure_exported"])
        self.assertTrue(all(self.payload["gatekeeper_checks"].values()))

    def test_docs(self) -> None:
        self.assertIn("stationary-refinement certificate", MD.read_text(encoding="utf-8"))
        self.assertIn("P2447/S1397", (ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md").read_text(encoding="utf-8"))
        self.assertIn("P2447/S1397", (ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md").read_text(encoding="utf-8"))


if __name__ == "__main__":
    unittest.main()
