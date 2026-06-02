import json
import subprocess
import sys
import unittest
from decimal import Decimal
from pathlib import Path

ROOT = Path(__file__).resolve().parent
SCRIPT = ROOT / "p2454_s1404_strict_pointwise_directed_decimal_weakest_band_replay_certificate.py"
OUT = ROOT / "generated" / "p2454_s1404_strict_pointwise_directed_decimal_weakest_band_replay_certificate.json"
MD = ROOT / "generated" / "p2454_s1404_strict_pointwise_directed_decimal_weakest_band_replay_certificate.md"
P2453 = ROOT / "generated" / "p2453_s1403_strict_pointwise_directed_decimal_weakest_cell_replay_certificate.json"


class P2454StrictPointwiseDirectedDecimalWeakestBandReplayCertificateTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        if not P2453.exists():
            subprocess.run([sys.executable, str(ROOT / "p2453_s1403_strict_pointwise_directed_decimal_weakest_cell_replay_certificate.py")], check=True)
        subprocess.run([sys.executable, str(SCRIPT)], check=True)
        cls.payload = json.loads(OUT.read_text(encoding="utf-8"))
        cls.cert = cls.payload["strict_pointwise_directed_decimal_weakest_band_replay_certificate"]
        cls.theorem = cls.cert["theorem_export"]

    def test_identity_and_band_counts(self) -> None:
        self.assertEqual(self.payload["packet_id"], "P2454")
        self.assertEqual(self.payload["stage_id"], "S1404")
        self.assertEqual(self.payload["status"], "PASS_STRICT_POINTWISE_DIRECTED_DECIMAL_WEAKEST_BAND_REPLAY_NO_FULL_INTERVAL_SELECTOR_SOURCE_THEOREM")
        self.assertTrue(OUT.exists())
        self.assertTrue(MD.exists())
        self.assertEqual(self.theorem["zero_projection_weakest_band_replay"]["band_cell_count"], 12)
        self.assertEqual(self.theorem["stationary_factor_weakest_band_replay"]["band_cell_count"], 12)

    def test_bands_exclude_zero(self) -> None:
        self.assertTrue(self.theorem["p2453_both_weakest_cells_exclude_zero"])
        self.assertTrue(self.theorem["zero_projection_weakest_band_replay"]["all_cells_exclude_zero"])
        self.assertTrue(self.theorem["stationary_factor_weakest_band_replay"]["all_cells_exclude_zero"])
        self.assertTrue(self.theorem["both_bands_exclude_zero_under_decimal_taylor_replay"])
        self.assertGreater(Decimal(self.theorem["minimum_decimal_separation_across_bands"]), Decimal("0"))

    def test_hard_limits(self) -> None:
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
        self.assertIn("directed-decimal weakest-band replay", MD.read_text(encoding="utf-8"))
        self.assertIn("P2454/S1404", (ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md").read_text(encoding="utf-8"))
        self.assertIn("P2454/S1404", (ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md").read_text(encoding="utf-8"))


if __name__ == "__main__":
    unittest.main()
