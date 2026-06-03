import json
import subprocess
import sys
import unittest
from decimal import Decimal
from pathlib import Path

ROOT = Path(__file__).resolve().parent
SCRIPT = ROOT / "p2466_s1416_strict_pointwise_interval_decimal_adaptive_descent_tail_boundary_ledger_certificate.py"
OUT = ROOT / "generated" / "p2466_s1416_strict_pointwise_interval_decimal_adaptive_descent_tail_boundary_ledger_certificate.json"
MD = ROOT / "generated" / "p2466_s1416_strict_pointwise_interval_decimal_adaptive_descent_tail_boundary_ledger_certificate.md"
P2465 = ROOT / "generated" / "p2465_s1415_strict_pointwise_interval_decimal_adaptive_descent_horizon_ledger_certificate.json"


class P2466StrictPointwiseIntervalDecimalAdaptiveDescentTailBoundaryLedgerCertificateTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        if not P2465.exists():
            subprocess.run([sys.executable, str(ROOT / "p2465_s1415_strict_pointwise_interval_decimal_adaptive_descent_horizon_ledger_certificate.py")], check=True)
        subprocess.run([sys.executable, str(SCRIPT)], check=True)
        cls.payload = json.loads(OUT.read_text(encoding="utf-8"))
        cls.cert = cls.payload["strict_pointwise_interval_decimal_adaptive_descent_tail_boundary_ledger_certificate"]
        cls.theorem = cls.cert["theorem_export"]

    def test_identity_and_tail_counts(self) -> None:
        self.assertEqual(self.payload["packet_id"], "P2466")
        self.assertEqual(self.payload["stage_id"], "S1416")
        self.assertEqual(
            self.payload["status"],
            "PASS_STRICT_POINTWISE_INTERVAL_DECIMAL_ADAPTIVE_DESCENT_TAIL_BOUNDARY_LEDGER_NO_FULL_COMPLEMENT_SELECTOR_SOURCE_THEOREM",
        )
        self.assertTrue(OUT.exists())
        self.assertTrue(MD.exists())
        self.assertEqual(self.theorem["p2465_total_descent_horizon_replay_count_inherited"], 64)
        self.assertEqual(len(self.theorem["family_descent_tail_boundary_replays"]), 2)
        self.assertEqual(self.theorem["total_tail_boundary_replay_count"], 6361)

    def test_tail_cells_exclude_zero_and_reach_boundary(self) -> None:
        self.assertTrue(self.theorem["p2465_all_descent_horizon_cells_exclude_zero_inherited"])
        self.assertTrue(self.theorem["all_tail_boundary_cells_exclude_zero"])
        self.assertGreater(Decimal(self.theorem["minimum_tail_boundary_decimal_separation"]), Decimal("0"))
        for family in self.theorem["family_descent_tail_boundary_replays"]:
            self.assertGreater(family["tail_boundary_replay_count"], 0)
            self.assertTrue(family["all_tail_boundary_cells_exclude_zero"])
            self.assertGreater(Decimal(family["minimum_tail_boundary_decimal_separation"]), Decimal("0"))
            self.assertTrue(all(row["decimal_interval_excludes_zero"] for row in family["tail_boundary_rows"]))
            if family["descent_direction"] < 0:
                self.assertEqual(family["boundary_endpoint_uncovered_index"], 0)
                self.assertEqual(family["boundary_side"], "left")
            else:
                self.assertEqual(
                    family["boundary_endpoint_uncovered_index"],
                    family["tail_boundary_rows"][-1]["uncovered_count"] - 1,
                )
                self.assertEqual(family["boundary_side"], "right")

    def test_tail_monotonic_diagnostic_is_consistent(self) -> None:
        for family in self.theorem["family_descent_tail_boundary_replays"]:
            previous = Decimal(family["inherited_p2465_endpoint_decimal_separation"])
            first_non_decreasing = None
            for row in family["tail_boundary_rows"]:
                current = Decimal(row["decimal_separation_from_zero"])
                if first_non_decreasing is None and current >= previous:
                    first_non_decreasing = row["tail_step"]
                previous = current
            self.assertEqual(family["first_non_decreasing_tail_step"], first_non_decreasing)
            self.assertEqual(family["local_bracket_found_in_tail"], first_non_decreasing is not None)
            self.assertEqual(
                family["strictly_decreasing_from_p2465_endpoint_to_boundary"],
                first_non_decreasing is None,
            )

    def test_hard_limits(self) -> None:
        self.assertTrue(self.theorem["adaptive_descent_tail_boundary_ledger_exported"])
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
        self.assertIn("descent tail-boundary ledger certificate", MD.read_text(encoding="utf-8"))
        self.assertIn("P2466/S1416", (ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md").read_text(encoding="utf-8"))
        self.assertIn("P2466/S1416", (ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md").read_text(encoding="utf-8"))


if __name__ == "__main__":
    unittest.main()
