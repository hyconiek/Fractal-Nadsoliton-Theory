import json
import subprocess
import sys
import unittest
from decimal import Decimal
from pathlib import Path

ROOT = Path(__file__).resolve().parent
SCRIPT = ROOT / "p2456_s1406_strict_pointwise_decimal_root_window_boundary_band_replay_certificate.py"
OUT = ROOT / "generated" / "p2456_s1406_strict_pointwise_decimal_root_window_boundary_band_replay_certificate.json"
MD = ROOT / "generated" / "p2456_s1406_strict_pointwise_decimal_root_window_boundary_band_replay_certificate.md"
P2455 = ROOT / "generated" / "p2455_s1405_strict_pointwise_directed_decimal_weakest_band_separation_monotonicity_certificate.json"


class P2456StrictPointwiseDecimalRootWindowBoundaryBandReplayCertificateTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        if not P2455.exists():
            subprocess.run([sys.executable, str(ROOT / "p2455_s1405_strict_pointwise_directed_decimal_weakest_band_separation_monotonicity_certificate.py")], check=True)
        subprocess.run([sys.executable, str(SCRIPT)], check=True)
        cls.payload = json.loads(OUT.read_text(encoding="utf-8"))
        cls.cert = cls.payload["strict_pointwise_decimal_root_window_boundary_band_replay_certificate"]
        cls.theorem = cls.cert["theorem_export"]

    def test_identity_and_boundary_counts(self) -> None:
        self.assertEqual(self.payload["packet_id"], "P2456")
        self.assertEqual(self.payload["stage_id"], "S1406")
        self.assertEqual(self.payload["status"], "PASS_STRICT_POINTWISE_DECIMAL_ROOT_WINDOW_BOUNDARY_BAND_REPLAY_NO_FULL_INTERVAL_SELECTOR_SOURCE_THEOREM")
        self.assertTrue(OUT.exists())
        self.assertTrue(MD.exists())
        self.assertTrue(self.theorem["p2455_bands_strictly_increasing_inherited"])
        self.assertEqual(self.theorem["total_boundary_band_count"], 6)
        self.assertEqual(self.theorem["total_replayed_boundary_cell_count"], 36)

    def test_boundary_bands_exclude_zero(self) -> None:
        self.assertEqual(len(self.theorem["zero_projection_boundary_band_replays"]), 4)
        self.assertEqual(len(self.theorem["stationary_factor_boundary_band_replays"]), 2)
        self.assertTrue(self.theorem["all_boundary_bands_exclude_zero"])
        self.assertGreater(Decimal(self.theorem["minimum_decimal_separation_across_all_boundary_bands"]), Decimal("0"))

    def test_hard_limits(self) -> None:
        self.assertTrue(self.theorem["bounded_all_root_window_boundary_replay_exported"])
        self.assertFalse(self.theorem["full_complement_directed_rounding_interval_theorem_exported_by_this_certificate"])
        self.assertFalse(self.theorem["symbolic_root_exclusion_theorem_exported_by_this_certificate"])
        self.assertFalse(self.theorem["pointwise_coordinate_selector_exported_by_this_certificate"])
        self.assertFalse(self.theorem["strict_observable_source_constraint_exported_by_this_certificate"])
        self.assertFalse(self.theorem["gauge_slice_theorem_exported_by_this_certificate"])
        self.assertFalse(self.theorem["strict_physical_value_generator_exported"])
        self.assertFalse(self.theorem["qw2191_discharged"])
        self.assertFalse(self.theorem["role_bearing_ltotal_exported"])
        self.assertFalse(self.theorem["toe_closure_exported"])
        self.assertTrue(all(self.payload["gatekeeper_checks"].values()))

    def test_docs(self) -> None:
        self.assertIn("root-window boundary-band replay certificate", MD.read_text(encoding="utf-8"))
        self.assertIn("P2456/S1406", (ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md").read_text(encoding="utf-8"))
        self.assertIn("P2456/S1406", (ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md").read_text(encoding="utf-8"))


if __name__ == "__main__":
    unittest.main()
