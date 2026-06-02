import json
import subprocess
import sys
import unittest
from decimal import Decimal
from pathlib import Path

ROOT = Path(__file__).resolve().parent
SCRIPT = ROOT / "p2460_s1410_strict_pointwise_interval_decimal_gap_sentinel_replay_certificate.py"
OUT = ROOT / "generated" / "p2460_s1410_strict_pointwise_interval_decimal_gap_sentinel_replay_certificate.json"
MD = ROOT / "generated" / "p2460_s1410_strict_pointwise_interval_decimal_gap_sentinel_replay_certificate.md"
P2459 = ROOT / "generated" / "p2459_s1409_strict_pointwise_interval_decimal_coverage_gap_ledger_certificate.json"


class P2460StrictPointwiseIntervalDecimalGapSentinelReplayCertificateTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        if not P2459.exists():
            subprocess.run([sys.executable, str(ROOT / "p2459_s1409_strict_pointwise_interval_decimal_coverage_gap_ledger_certificate.py")], check=True)
        subprocess.run([sys.executable, str(SCRIPT)], check=True)
        cls.payload = json.loads(OUT.read_text(encoding="utf-8"))
        cls.cert = cls.payload["strict_pointwise_interval_decimal_gap_sentinel_replay_certificate"]
        cls.theorem = cls.cert["theorem_export"]

    def test_identity_and_gap_counts(self) -> None:
        self.assertEqual(self.payload["packet_id"], "P2460")
        self.assertEqual(self.payload["stage_id"], "S1410")
        self.assertEqual(
            self.payload["status"],
            "PASS_STRICT_POINTWISE_INTERVAL_DECIMAL_GAP_SENTINEL_REPLAY_NO_FULL_COMPLEMENT_SELECTOR_SOURCE_THEOREM",
        )
        self.assertTrue(OUT.exists())
        self.assertTrue(MD.exists())
        self.assertEqual(self.theorem["p2459_total_unreplayed_gap_inherited"], 99846)
        self.assertEqual(self.theorem["total_reconstructed_uncovered_gap_cell_count"], 99846)
        self.assertEqual(self.theorem["total_gap_sentinel_replay_count"], 25)

    def test_gap_sentinels_exclude_zero(self) -> None:
        self.assertTrue(self.theorem["all_gap_sentinels_exclude_zero"])
        self.assertGreater(Decimal(self.theorem["minimum_gap_sentinel_decimal_separation"]), Decimal("0"))
        self.assertEqual(len(self.theorem["family_gap_sentinel_replays"]), 2)
        for family in self.theorem["family_gap_sentinel_replays"]:
            self.assertTrue(family["all_gap_sentinels_exclude_zero"])
            self.assertGreater(family["uncovered_by_boundary_chain_cell_count"], 0)
            for segment in family["segment_rows"]:
                self.assertTrue(segment["all_segment_sentinels_exclude_zero"])
                self.assertGreater(segment["uncovered_by_boundary_chain_cell_count"], 0)
                self.assertGreaterEqual(segment["sentinel_replay_count"], 5)

    def test_hard_limits(self) -> None:
        self.assertTrue(self.theorem["gap_sentinel_replay_exported"])
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
        self.assertIn("gap-sentinel replay certificate", MD.read_text(encoding="utf-8"))
        self.assertIn("P2460/S1410", (ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md").read_text(encoding="utf-8"))
        self.assertIn("P2460/S1410", (ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md").read_text(encoding="utf-8"))


if __name__ == "__main__":
    unittest.main()
