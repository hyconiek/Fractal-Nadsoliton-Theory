import json
import subprocess
import sys
import unittest
from decimal import Decimal
from pathlib import Path

ROOT = Path(__file__).resolve().parent
SCRIPT = ROOT / "p2465_s1415_strict_pointwise_interval_decimal_adaptive_descent_horizon_ledger_certificate.py"
OUT = ROOT / "generated" / "p2465_s1415_strict_pointwise_interval_decimal_adaptive_descent_horizon_ledger_certificate.json"
MD = ROOT / "generated" / "p2465_s1415_strict_pointwise_interval_decimal_adaptive_descent_horizon_ledger_certificate.md"
P2464 = ROOT / "generated" / "p2464_s1414_strict_pointwise_interval_decimal_adaptive_flank_descent_extension_replay_certificate.json"


class P2465StrictPointwiseIntervalDecimalAdaptiveDescentHorizonLedgerCertificateTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        if not P2464.exists():
            subprocess.run([sys.executable, str(ROOT / "p2464_s1414_strict_pointwise_interval_decimal_adaptive_flank_descent_extension_replay_certificate.py")], check=True)
        subprocess.run([sys.executable, str(SCRIPT)], check=True)
        cls.payload = json.loads(OUT.read_text(encoding="utf-8"))
        cls.cert = cls.payload["strict_pointwise_interval_decimal_adaptive_descent_horizon_ledger_certificate"]
        cls.theorem = cls.cert["theorem_export"]

    def test_identity_and_horizon_counts(self) -> None:
        self.assertEqual(self.payload["packet_id"], "P2465")
        self.assertEqual(self.payload["stage_id"], "S1415")
        self.assertEqual(
            self.payload["status"],
            "PASS_STRICT_POINTWISE_INTERVAL_DECIMAL_ADAPTIVE_DESCENT_HORIZON_LEDGER_NO_FULL_COMPLEMENT_SELECTOR_SOURCE_THEOREM",
        )
        self.assertTrue(OUT.exists())
        self.assertTrue(MD.exists())
        self.assertEqual(self.theorem["p2464_total_descent_extension_replay_count_inherited"], 8)
        self.assertEqual(self.theorem["descent_horizon_steps_requested"], 32)
        self.assertEqual(len(self.theorem["family_descent_horizon_replays"]), 2)
        self.assertEqual(self.theorem["total_descent_horizon_replay_count"], 64)

    def test_horizon_cells_exclude_zero_and_remain_unbracketed(self) -> None:
        self.assertTrue(self.theorem["p2464_all_descent_extension_cells_exclude_zero_inherited"])
        self.assertTrue(self.theorem["all_descent_horizon_cells_exclude_zero"])
        self.assertTrue(self.theorem["all_horizons_strictly_decreasing_from_p2464_endpoint"])
        self.assertFalse(self.theorem["any_local_bracket_found_within_horizon"])
        self.assertGreater(Decimal(self.theorem["minimum_descent_horizon_decimal_separation"]), Decimal("0"))
        for family in self.theorem["family_descent_horizon_replays"]:
            self.assertEqual(family["descent_horizon_replay_count"], 32)
            self.assertTrue(family["all_descent_horizon_cells_exclude_zero"])
            self.assertTrue(family["strictly_decreasing_from_p2464_endpoint_along_horizon"])
            self.assertFalse(family["local_bracket_found_within_horizon"])
            self.assertIsNone(family["first_non_decreasing_horizon_step"])
            previous = Decimal(family["inherited_p2464_endpoint_decimal_separation"])
            for row in family["descent_horizon_rows"]:
                current = Decimal(row["decimal_separation_from_zero"])
                self.assertTrue(row["decimal_interval_excludes_zero"])
                self.assertLess(current, previous)
                previous = current

    def test_hard_limits(self) -> None:
        self.assertTrue(self.theorem["adaptive_descent_horizon_ledger_exported"])
        self.assertFalse(self.theorem["decimal_full_complement_replay_exported_by_this_certificate"])
        self.assertFalse(self.theorem["directed_rounding_interval_theorem_exported_by_this_certificate"])
        self.assertFalse(self.theorem["symbolic_root_exclusion_theorem_exported_by_this_certificate"])
        self.assertFalse(self.theorem["analytic_monotonicity_theorem_exported_by_this_certificate"])
        self.assertFalse(self.theorem["pointwise_coordinate_selector_exported_by_this_certificate"])
        self.assertFalse(self.theorem["strict_observable_source_constraint_exported_by_this_certificate"])
        self.assertFalse(self.theorem["gauge_slice_theorem_exported_by_this_certificate"])
        self.assertFalse(self.theorem["strict_physical_value_generator_exported"])
        self.assertFalse(self.theorem["qw2191_discharged"])
        self.assertFalse(self.theorem["role_bearing_ltotal_exported"])
        self.assertFalse(self.theorem["toe_closure_exported"])
        self.assertTrue(all(self.payload["gatekeeper_checks"].values()))

    def test_docs(self) -> None:
        self.assertIn("adaptive descent-horizon ledger certificate", MD.read_text(encoding="utf-8"))
        self.assertIn("P2465/S1415", (ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md").read_text(encoding="utf-8"))
        self.assertIn("P2465/S1415", (ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md").read_text(encoding="utf-8"))


if __name__ == "__main__":
    unittest.main()
