import json
import subprocess
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent
SCRIPT = ROOT / "p2450_s1400_strict_pointwise_projection_root_isolation_margin_certificate.py"
OUT = ROOT / "generated" / "p2450_s1400_strict_pointwise_projection_root_isolation_margin_certificate.json"
MD = ROOT / "generated" / "p2450_s1400_strict_pointwise_projection_root_isolation_margin_certificate.md"
P2449 = ROOT / "generated" / "p2449_s1399_strict_pointwise_rank_lift_projection_reduction_certificate.json"


class P2450StrictPointwiseProjectionRootIsolationMarginCertificateTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        if not P2449.exists():
            subprocess.run([sys.executable, str(ROOT / "p2449_s1399_strict_pointwise_rank_lift_projection_reduction_certificate.py")], check=True)
        subprocess.run([sys.executable, str(SCRIPT)], check=True)
        cls.payload = json.loads(OUT.read_text(encoding="utf-8"))
        cls.cert = cls.payload["strict_pointwise_projection_root_isolation_margin_certificate"]
        cls.theorem = cls.cert["theorem_export"]
        cls.zero_projection = cls.theorem["zero_projection_amplitude_certificate"]
        cls.stationary_factor = cls.theorem["stationary_factor_certificate"]

    def test_identity_outputs_and_inheritance(self) -> None:
        self.assertEqual(self.payload["packet_id"], "P2450")
        self.assertEqual(self.payload["stage_id"], "S1400")
        self.assertEqual(self.payload["status"], "PASS_STRICT_POINTWISE_PROJECTION_ROOT_ISOLATION_MARGIN_NO_EXACT_INTERVAL_SELECTOR_SOURCE_THEOREM")
        self.assertTrue(OUT.exists())
        self.assertTrue(MD.exists())
        self.assertEqual(self.zero_projection["root_count"], 2)
        self.assertEqual(self.stationary_factor["root_count"], 1)

    def test_sampled_isolation_margins(self) -> None:
        self.assertTrue(self.theorem["all_root_windows_sign_change"])
        self.assertTrue(self.theorem["all_root_windows_sampled_monotone"])
        self.assertTrue(self.theorem["all_complement_exclusion_margins_positive"])
        self.assertGreater(self.theorem["minimum_exclusion_margin_across_families"], 0.0)
        self.assertEqual(self.zero_projection["sampled_lipschitz_exclusion"]["failed_cell_count"], 0)
        self.assertEqual(self.stationary_factor["sampled_lipschitz_exclusion"]["failed_cell_count"], 0)
        self.assertGreater(self.zero_projection["sampled_lipschitz_exclusion"]["cell_count"], 4000)
        self.assertGreater(self.stationary_factor["sampled_lipschitz_exclusion"]["cell_count"], 4000)
        self.assertTrue(self.payload["gatekeeper_checks"]["best_root_volume_replayed"])

    def test_hard_limits(self) -> None:
        self.assertTrue(self.theorem["sampled_lipschitz_margin_certificate_exported"])
        self.assertFalse(self.theorem["exact_interval_root_exclusion_theorem_exported_by_this_certificate"])
        self.assertFalse(self.theorem["pointwise_coordinate_selector_exported_by_this_certificate"])
        self.assertFalse(self.theorem["strict_observable_source_constraint_exported_by_this_certificate"])
        self.assertFalse(self.theorem["gauge_slice_theorem_exported_by_this_certificate"])
        self.assertFalse(self.theorem["strict_physical_value_generator_exported"])
        self.assertFalse(self.theorem["qw2191_discharged"])
        self.assertFalse(self.theorem["role_bearing_ltotal_exported"])
        self.assertFalse(self.theorem["toe_closure_exported"])
        self.assertTrue(all(self.payload["gatekeeper_checks"].values()))

    def test_docs(self) -> None:
        self.assertIn("root-isolation margin certificate", MD.read_text(encoding="utf-8"))
        self.assertIn("P2450/S1400", (ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md").read_text(encoding="utf-8"))
        self.assertIn("P2450/S1400", (ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md").read_text(encoding="utf-8"))


if __name__ == "__main__":
    unittest.main()
