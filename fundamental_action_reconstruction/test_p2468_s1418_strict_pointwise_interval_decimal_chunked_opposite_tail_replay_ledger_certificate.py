import json
import subprocess
import sys
import unittest
from decimal import Decimal
from pathlib import Path

ROOT = Path(__file__).resolve().parent
SCRIPT = ROOT / "p2468_s1418_strict_pointwise_interval_decimal_chunked_opposite_tail_replay_ledger_certificate.py"
OUT = ROOT / "generated" / "p2468_s1418_strict_pointwise_interval_decimal_chunked_opposite_tail_replay_ledger_certificate.json"
MD = ROOT / "generated" / "p2468_s1418_strict_pointwise_interval_decimal_chunked_opposite_tail_replay_ledger_certificate.md"
P2467 = ROOT / "generated" / "p2467_s1417_strict_pointwise_interval_decimal_opposite_tail_sentinel_ledger_certificate.json"


class P2468StrictPointwiseIntervalDecimalChunkedOppositeTailReplayLedgerCertificateTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        if not P2467.exists():
            subprocess.run([sys.executable, str(ROOT / "p2467_s1417_strict_pointwise_interval_decimal_opposite_tail_sentinel_ledger_certificate.py")], check=True)
        subprocess.run([sys.executable, str(SCRIPT)], check=True)
        cls.payload = json.loads(OUT.read_text(encoding="utf-8"))
        cls.cert = cls.payload["strict_pointwise_interval_decimal_chunked_opposite_tail_replay_ledger_certificate"]
        cls.theorem = cls.cert["theorem_export"]

    def test_identity_and_inherited_counts(self) -> None:
        self.assertEqual(self.payload["packet_id"], "P2468")
        self.assertEqual(self.payload["stage_id"], "S1418")
        self.assertEqual(
            self.payload["status"],
            "PASS_STRICT_POINTWISE_INTERVAL_DECIMAL_CHUNKED_OPPOSITE_TAIL_REPLAY_LEDGER_NO_FULL_TAIL_COMPLEMENT_SELECTOR_SOURCE_THEOREM",
        )
        self.assertTrue(OUT.exists())
        self.assertTrue(MD.exists())
        self.assertEqual(self.theorem["p2451_total_interval_complement_cell_count_inherited"], 99882)
        self.assertEqual(self.theorem["p2459_total_unreplayed_by_decimal_boundary_chain_cell_count_inherited"], 99846)
        self.assertEqual(self.theorem["p2466_total_tail_boundary_replay_count_inherited"], 6361)
        self.assertEqual(self.theorem["p2467_total_opposite_tail_available_cell_count_inherited"], 45165)
        self.assertEqual(self.theorem["p2467_total_opposite_tail_sentinel_replay_count_inherited"], 42)

    def test_chunked_replay_counts_and_budget(self) -> None:
        self.assertEqual(len(self.theorem["family_chunked_opposite_tail_replays"]), 2)
        self.assertEqual(self.theorem["total_opposite_tail_available_cell_count"], 45165)
        self.assertEqual(self.theorem["total_opposite_tail_chunk_count"], 46)
        self.assertEqual(self.theorem["total_opposite_tail_chunked_replay_count"], 140)
        self.assertEqual(self.theorem["total_opposite_tail_unreplayed_budget_after_chunked_replay"], 45025)
        self.assertEqual(self.theorem["p2459_unreplayed_budget_not_covered_by_p2466_p2468_tail_replays"], 93345)
        for family in self.theorem["family_chunked_opposite_tail_replays"]:
            self.assertGreater(family["opposite_tail_available_count"], family["opposite_tail_chunked_replay_count"])
            self.assertGreater(family["opposite_tail_unreplayed_budget_after_chunked_replay"], 0)
            self.assertFalse(family["full_opposite_tail_replay_exported_by_this_family"])
            self.assertEqual(
                family["opposite_tail_available_count"] - family["opposite_tail_chunked_replay_count"],
                family["opposite_tail_unreplayed_budget_after_chunked_replay"],
            )

    def test_replayed_cells_exclude_zero_and_are_disjoint(self) -> None:
        self.assertTrue(self.theorem["p2466_all_tail_boundary_cells_exclude_zero_inherited"])
        self.assertTrue(self.theorem["p2467_all_opposite_tail_sentinels_exclude_zero_inherited"])
        self.assertTrue(self.theorem["all_opposite_tail_chunked_replayed_cells_exclude_zero"])
        self.assertTrue(self.theorem["all_opposite_tail_chunked_replay_indexes_unique_by_family"])
        self.assertTrue(self.theorem["all_opposite_tail_chunked_replays_disjoint_from_p2466_descent_tail"])
        self.assertGreater(Decimal(self.theorem["minimum_opposite_tail_chunked_replay_decimal_separation"]), Decimal("0"))
        for family in self.theorem["family_chunked_opposite_tail_replays"]:
            self.assertTrue(family["all_opposite_tail_chunked_replayed_cells_exclude_zero"])
            self.assertTrue(family["all_opposite_tail_chunked_replay_indexes_unique"])
            self.assertTrue(family["all_opposite_tail_chunked_replays_disjoint_from_p2466_descent_tail"])
            self.assertGreater(Decimal(family["minimum_opposite_tail_chunked_replay_decimal_separation"]), Decimal("0"))
            self.assertTrue(all(chunk["all_chunk_replayed_cells_exclude_zero"] for chunk in family["opposite_tail_chunk_replays"]))

    def test_hard_limits_and_docs(self) -> None:
        self.assertTrue(self.theorem["chunked_opposite_tail_replay_ledger_exported"])
        self.assertFalse(self.theorem["full_opposite_tail_replay_exported_by_this_certificate"])
        self.assertFalse(self.theorem["remaining_complement_segments_replayed_by_this_certificate"])
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
        self.assertIn("chunked opposite-tail replay ledger", MD.read_text(encoding="utf-8"))
        self.assertIn("P2468/S1418", (ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md").read_text(encoding="utf-8"))
        self.assertIn("P2468/S1418", (ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md").read_text(encoding="utf-8"))


if __name__ == "__main__":
    unittest.main()
