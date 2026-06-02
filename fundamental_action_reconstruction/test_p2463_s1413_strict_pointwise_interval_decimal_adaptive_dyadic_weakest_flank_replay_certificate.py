import json
import subprocess
import sys
import unittest
from decimal import Decimal
from pathlib import Path

ROOT = Path(__file__).resolve().parent
SCRIPT = ROOT / "p2463_s1413_strict_pointwise_interval_decimal_adaptive_dyadic_weakest_flank_replay_certificate.py"
OUT = ROOT / "generated" / "p2463_s1413_strict_pointwise_interval_decimal_adaptive_dyadic_weakest_flank_replay_certificate.json"
MD = ROOT / "generated" / "p2463_s1413_strict_pointwise_interval_decimal_adaptive_dyadic_weakest_flank_replay_certificate.md"
P2462 = ROOT / "generated" / "p2462_s1412_strict_pointwise_interval_decimal_gap_dyadic_refinement_replay_certificate.json"


class P2463StrictPointwiseIntervalDecimalAdaptiveDyadicWeakestFlankReplayCertificateTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        if not P2462.exists():
            subprocess.run([sys.executable, str(ROOT / "p2462_s1412_strict_pointwise_interval_decimal_gap_dyadic_refinement_replay_certificate.py")], check=True)
        subprocess.run([sys.executable, str(SCRIPT)], check=True)
        cls.payload = json.loads(OUT.read_text(encoding="utf-8"))
        cls.cert = cls.payload["strict_pointwise_interval_decimal_adaptive_dyadic_weakest_flank_replay_certificate"]
        cls.theorem = cls.cert["theorem_export"]

    def test_identity_and_adaptive_counts(self) -> None:
        self.assertEqual(self.payload["packet_id"], "P2463")
        self.assertEqual(self.payload["stage_id"], "S1413")
        self.assertEqual(
            self.payload["status"],
            "PASS_STRICT_POINTWISE_INTERVAL_DECIMAL_ADAPTIVE_DYADIC_WEAKEST_FLANK_REPLAY_NO_FULL_COMPLEMENT_SELECTOR_SOURCE_THEOREM",
        )
        self.assertTrue(OUT.exists())
        self.assertTrue(MD.exists())
        self.assertEqual(self.theorem["p2460_total_gap_sentinel_replay_count_inherited"], 25)
        self.assertEqual(self.theorem["p2462_total_dyadic_refinement_replay_count_inherited"], 20)
        self.assertEqual(self.theorem["adaptive_flank_radius_requested"], 4)
        self.assertEqual(len(self.theorem["family_adaptive_weakest_flank_replays"]), 2)
        self.assertEqual(self.theorem["total_adaptive_flank_replay_count"], 16)

    def test_adaptive_flanks_exclude_zero_and_find_descent(self) -> None:
        self.assertTrue(self.theorem["p2462_all_dyadic_refinement_cells_exclude_zero_inherited"])
        self.assertTrue(self.theorem["all_adaptive_flank_cells_exclude_zero"])
        self.assertTrue(self.theorem["all_families_found_smaller_than_inherited_weakest_dyadic"])
        self.assertGreater(self.theorem["total_smaller_than_inherited_weakest_dyadic_count"], 0)
        self.assertGreater(Decimal(self.theorem["minimum_adaptive_flank_decimal_separation"]), Decimal("0"))
        for family in self.theorem["family_adaptive_weakest_flank_replays"]:
            self.assertEqual(family["adaptive_flank_replay_count"], 8)
            self.assertTrue(family["all_adaptive_flank_cells_exclude_zero"])
            self.assertTrue(family["found_smaller_than_inherited_weakest_dyadic"])
            self.assertGreater(family["smaller_than_inherited_weakest_dyadic_count"], 0)
            self.assertLess(
                Decimal(family["minimum_adaptive_flank_decimal_separation"]),
                Decimal(family["inherited_weakest_dyadic_decimal_separation"]),
            )
            for row in family["adaptive_flank_rows"]:
                self.assertTrue(row["decimal_interval_excludes_zero"])
                self.assertNotEqual(row["offset_from_p2462_weakest_dyadic"], 0)

    def test_hard_limits(self) -> None:
        self.assertTrue(self.theorem["adaptive_dyadic_weakest_flank_replay_exported"])
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
        self.assertIn("adaptive dyadic weakest-flank replay certificate", MD.read_text(encoding="utf-8"))
        self.assertIn("P2463/S1413", (ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md").read_text(encoding="utf-8"))
        self.assertIn("P2463/S1413", (ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md").read_text(encoding="utf-8"))


if __name__ == "__main__":
    unittest.main()
