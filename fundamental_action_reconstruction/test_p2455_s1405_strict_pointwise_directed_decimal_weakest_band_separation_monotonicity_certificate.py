import json
import subprocess
import sys
import unittest
from decimal import Decimal
from pathlib import Path

ROOT = Path(__file__).resolve().parent
SCRIPT = ROOT / "p2455_s1405_strict_pointwise_directed_decimal_weakest_band_separation_monotonicity_certificate.py"
OUT = ROOT / "generated" / "p2455_s1405_strict_pointwise_directed_decimal_weakest_band_separation_monotonicity_certificate.json"
MD = ROOT / "generated" / "p2455_s1405_strict_pointwise_directed_decimal_weakest_band_separation_monotonicity_certificate.md"
P2454 = ROOT / "generated" / "p2454_s1404_strict_pointwise_directed_decimal_weakest_band_replay_certificate.json"


class P2455StrictPointwiseDirectedDecimalWeakestBandSeparationMonotonicityCertificateTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        if not P2454.exists():
            subprocess.run([sys.executable, str(ROOT / "p2454_s1404_strict_pointwise_directed_decimal_weakest_band_replay_certificate.py")], check=True)
        subprocess.run([sys.executable, str(SCRIPT)], check=True)
        cls.payload = json.loads(OUT.read_text(encoding="utf-8"))
        cls.cert = cls.payload["strict_pointwise_directed_decimal_weakest_band_separation_monotonicity_certificate"]
        cls.theorem = cls.cert["theorem_export"]

    def test_identity_and_inheritance(self) -> None:
        self.assertEqual(self.payload["packet_id"], "P2455")
        self.assertEqual(self.payload["stage_id"], "S1405")
        self.assertEqual(self.payload["status"], "PASS_STRICT_POINTWISE_DIRECTED_DECIMAL_WEAKEST_BAND_SEPARATION_MONOTONICITY_NO_FULL_INTERVAL_SELECTOR_SOURCE_THEOREM")
        self.assertTrue(OUT.exists())
        self.assertTrue(MD.exists())
        self.assertTrue(self.theorem["inherited_p2454_bands_exclude_zero"])

    def test_separation_monotonicity(self) -> None:
        zero = self.theorem["zero_projection_separation_monotonicity_audit"]
        stationary = self.theorem["stationary_factor_separation_monotonicity_audit"]
        self.assertTrue(zero["strictly_increasing_separation"])
        self.assertTrue(stationary["strictly_increasing_separation"])
        self.assertEqual(zero["weakest_decimal_separation_index"], 0)
        self.assertEqual(stationary["weakest_decimal_separation_index"], 0)
        self.assertTrue(self.theorem["both_bands_have_strictly_increasing_separation"])
        self.assertTrue(self.theorem["both_bands_weakest_at_boundary_cell"])
        self.assertGreater(Decimal(self.theorem["minimum_first_difference_across_bands"]), Decimal("0"))

    def test_hard_limits(self) -> None:
        self.assertTrue(self.theorem["bounded_critical_band_monotonicity_certificate_exported"])
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
        self.assertIn("separation-monotonicity certificate", MD.read_text(encoding="utf-8"))
        self.assertIn("P2455/S1405", (ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md").read_text(encoding="utf-8"))
        self.assertIn("P2455/S1405", (ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md").read_text(encoding="utf-8"))


if __name__ == "__main__":
    unittest.main()
