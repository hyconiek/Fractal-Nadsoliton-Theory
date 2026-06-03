import json
import subprocess
import sys
import unittest
from decimal import Decimal
from pathlib import Path

ROOT = Path(__file__).resolve().parent
SCRIPT = ROOT / "p2467_s1417_strict_pointwise_interval_decimal_opposite_tail_sentinel_ledger_certificate.py"
OUT = ROOT / "generated" / "p2467_s1417_strict_pointwise_interval_decimal_opposite_tail_sentinel_ledger_certificate.json"
MD = ROOT / "generated" / "p2467_s1417_strict_pointwise_interval_decimal_opposite_tail_sentinel_ledger_certificate.md"
P2466 = ROOT / "generated" / "p2466_s1416_strict_pointwise_interval_decimal_adaptive_descent_tail_boundary_ledger_certificate.json"


class P2467StrictPointwiseIntervalDecimalOppositeTailSentinelLedgerCertificateTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        if not P2466.exists():
            subprocess.run([sys.executable, str(ROOT / "p2466_s1416_strict_pointwise_interval_decimal_adaptive_descent_tail_boundary_ledger_certificate.py")], check=True)
        subprocess.run([sys.executable, str(SCRIPT)], check=True)
        cls.payload = json.loads(OUT.read_text(encoding="utf-8"))
        cls.cert = cls.payload["strict_pointwise_interval_decimal_opposite_tail_sentinel_ledger_certificate"]
        cls.theorem = cls.cert["theorem_export"]

    def test_identity_and_sentinel_counts(self) -> None:
        self.assertEqual(self.payload["packet_id"], "P2467")
        self.assertEqual(self.payload["stage_id"], "S1417")
        self.assertEqual(
            self.payload["status"],
            "PASS_STRICT_POINTWISE_INTERVAL_DECIMAL_OPPOSITE_TAIL_SENTINEL_LEDGER_NO_FULL_TAIL_COMPLEMENT_SELECTOR_SOURCE_THEOREM",
        )
        self.assertTrue(OUT.exists())
        self.assertTrue(MD.exists())
        self.assertEqual(self.theorem["p2466_total_tail_boundary_replay_count_inherited"], 6361)
        self.assertEqual(len(self.theorem["family_opposite_tail_sentinel_replays"]), 2)
        self.assertEqual(self.theorem["total_opposite_tail_available_cell_count"], 45165)
        self.assertEqual(self.theorem["total_opposite_tail_sentinel_replay_count"], 42)

    def test_sentinels_exclude_zero_and_are_disjoint(self) -> None:
        self.assertTrue(self.theorem["p2466_all_tail_boundary_cells_exclude_zero_inherited"])
        self.assertTrue(self.theorem["all_opposite_tail_sentinels_exclude_zero"])
        self.assertTrue(self.theorem["all_opposite_tail_sentinel_indexes_unique"])
        self.assertTrue(self.theorem["all_opposite_tail_sentinels_disjoint_from_p2466_descent_tail"])
        self.assertGreater(Decimal(self.theorem["minimum_opposite_tail_sentinel_decimal_separation"]), Decimal("0"))
        for family in self.theorem["family_opposite_tail_sentinel_replays"]:
            self.assertGreater(family["opposite_tail_available_count"], family["opposite_tail_sentinel_replay_count"])
            self.assertEqual(family["opposite_tail_sentinel_replay_count"], 21)
            self.assertTrue(family["all_opposite_tail_sentinels_exclude_zero"])
            self.assertTrue(family["all_opposite_tail_sentinel_indexes_unique"])
            self.assertTrue(family["all_opposite_tail_sentinels_disjoint_from_p2466_descent_tail"])
            self.assertGreater(Decimal(family["minimum_opposite_tail_sentinel_decimal_separation"]), Decimal("0"))
            indexes = [row["uncovered_index"] for row in family["opposite_tail_sentinel_rows"]]
            self.assertEqual(len(indexes), len(set(indexes)))
            self.assertTrue(all(row["decimal_interval_excludes_zero"] for row in family["opposite_tail_sentinel_rows"]))
            self.assertEqual(family["opposite_direction"], -family["inherited_p2466_descent_direction"])

    def test_hard_limits(self) -> None:
        self.assertTrue(self.theorem["opposite_tail_sentinel_ledger_exported"])
        self.assertFalse(self.theorem["opposite_tail_full_replay_exported_by_this_certificate"])
        self.assertFalse(self.theorem["remaining_complement_segments_replayed_by_this_certificate"])
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
        self.assertIn("opposite-tail sentinel ledger certificate", MD.read_text(encoding="utf-8"))
        self.assertIn("P2467/S1417", (ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md").read_text(encoding="utf-8"))
        self.assertIn("P2467/S1417", (ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md").read_text(encoding="utf-8"))


if __name__ == "__main__":
    unittest.main()
