import json
import subprocess
import sys
import unittest
from decimal import Decimal
from pathlib import Path

ROOT = Path(__file__).resolve().parent
SCRIPT = ROOT / "p2462_s1412_strict_pointwise_interval_decimal_gap_dyadic_refinement_replay_certificate.py"
OUT = ROOT / "generated" / "p2462_s1412_strict_pointwise_interval_decimal_gap_dyadic_refinement_replay_certificate.json"
MD = ROOT / "generated" / "p2462_s1412_strict_pointwise_interval_decimal_gap_dyadic_refinement_replay_certificate.md"
P2460 = ROOT / "generated" / "p2460_s1410_strict_pointwise_interval_decimal_gap_sentinel_replay_certificate.json"
P2461 = ROOT / "generated" / "p2461_s1411_strict_pointwise_interval_decimal_gap_weakest_neighborhood_replay_certificate.json"


class P2462StrictPointwiseIntervalDecimalGapDyadicRefinementReplayCertificateTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        if not P2460.exists():
            subprocess.run([sys.executable, str(ROOT / "p2460_s1410_strict_pointwise_interval_decimal_gap_sentinel_replay_certificate.py")], check=True)
        if not P2461.exists():
            subprocess.run([sys.executable, str(ROOT / "p2461_s1411_strict_pointwise_interval_decimal_gap_weakest_neighborhood_replay_certificate.py")], check=True)
        subprocess.run([sys.executable, str(SCRIPT)], check=True)
        cls.payload = json.loads(OUT.read_text(encoding="utf-8"))
        cls.cert = cls.payload["strict_pointwise_interval_decimal_gap_dyadic_refinement_replay_certificate"]
        cls.theorem = cls.cert["theorem_export"]

    def test_identity_and_dyadic_counts(self) -> None:
        self.assertEqual(self.payload["packet_id"], "P2462")
        self.assertEqual(self.payload["stage_id"], "S1412")
        self.assertEqual(
            self.payload["status"],
            "PASS_STRICT_POINTWISE_INTERVAL_DECIMAL_GAP_DYADIC_REFINEMENT_REPLAY_NO_FULL_COMPLEMENT_SELECTOR_SOURCE_THEOREM",
        )
        self.assertTrue(OUT.exists())
        self.assertTrue(MD.exists())
        self.assertEqual(self.theorem["p2460_total_gap_sentinel_replay_count_inherited"], 25)
        self.assertEqual(self.theorem["p2461_total_neighborhood_replay_count_inherited"], 8)
        self.assertEqual(self.theorem["dyadic_refinement_fractions"], ["0.125", "0.375", "0.625", "0.875"])
        self.assertEqual(len(self.theorem["family_dyadic_refinement_replays"]), 2)
        self.assertEqual(self.theorem["total_dyadic_refinement_replay_count"], 20)

    def test_dyadic_refinements_exclude_zero_and_avoid_quarter_sentinels(self) -> None:
        self.assertTrue(self.theorem["p2460_all_gap_sentinels_exclude_zero_inherited"])
        self.assertTrue(self.theorem["p2461_all_neighborhood_cells_exclude_zero_inherited"])
        self.assertTrue(self.theorem["all_dyadic_refinement_cells_exclude_zero"])
        self.assertTrue(self.theorem["all_refinement_indexes_avoid_p2460_quarter_sentinels"])
        self.assertGreater(Decimal(self.theorem["minimum_dyadic_refinement_decimal_separation"]), Decimal("0"))
        for family in self.theorem["family_dyadic_refinement_replays"]:
            self.assertEqual(family["dyadic_refinement_replay_count"], 4 * family["segment_count"])
            self.assertTrue(family["all_dyadic_refinement_cells_exclude_zero"])
            self.assertTrue(family["all_refinement_indexes_avoid_p2460_quarter_sentinels"])
            for segment in family["segment_rows"]:
                self.assertEqual(segment["dyadic_refinement_replay_count"], 4)
                self.assertTrue(segment["no_p2460_quarter_sentinel_index_duplicates"])
                for row in segment["dyadic_refinement_replays"]:
                    self.assertTrue(row["decimal_interval_excludes_zero"])
                    self.assertIn(row["refinement_fraction"], self.theorem["dyadic_refinement_fractions"])

    def test_hard_limits(self) -> None:
        self.assertTrue(self.theorem["gap_dyadic_refinement_replay_exported"])
        self.assertFalse(self.theorem["decimal_full_complement_replay_exported_by_this_certificate"])
        self.assertFalse(self.theorem["directed_rounding_interval_theorem_exported_by_this_certificate"])
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
        self.assertIn("gap dyadic-refinement replay certificate", MD.read_text(encoding="utf-8"))
        self.assertIn("P2462/S1412", (ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md").read_text(encoding="utf-8"))
        self.assertIn("P2462/S1412", (ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md").read_text(encoding="utf-8"))


if __name__ == "__main__":
    unittest.main()
