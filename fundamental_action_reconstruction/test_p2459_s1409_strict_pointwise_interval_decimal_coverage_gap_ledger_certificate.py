import json
import subprocess
import sys
import unittest
from decimal import Decimal
from pathlib import Path

ROOT = Path(__file__).resolve().parent
SCRIPT = ROOT / "p2459_s1409_strict_pointwise_interval_decimal_coverage_gap_ledger_certificate.py"
OUT = ROOT / "generated" / "p2459_s1409_strict_pointwise_interval_decimal_coverage_gap_ledger_certificate.json"
MD = ROOT / "generated" / "p2459_s1409_strict_pointwise_interval_decimal_coverage_gap_ledger_certificate.md"
P2458 = ROOT / "generated" / "p2458_s1408_strict_pointwise_interval_decimal_weakest_cell_alignment_certificate.json"


class P2459StrictPointwiseIntervalDecimalCoverageGapLedgerCertificateTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        if not P2458.exists():
            subprocess.run([sys.executable, str(ROOT / "p2458_s1408_strict_pointwise_interval_decimal_weakest_cell_alignment_certificate.py")], check=True)
        subprocess.run([sys.executable, str(SCRIPT)], check=True)
        cls.payload = json.loads(OUT.read_text(encoding="utf-8"))
        cls.cert = cls.payload["strict_pointwise_interval_decimal_coverage_gap_ledger_certificate"]
        cls.theorem = cls.cert["theorem_export"]

    def test_identity_and_coverage_counts(self) -> None:
        self.assertEqual(self.payload["packet_id"], "P2459")
        self.assertEqual(self.payload["stage_id"], "S1409")
        self.assertEqual(
            self.payload["status"],
            "PASS_STRICT_POINTWISE_INTERVAL_DECIMAL_COVERAGE_GAP_LEDGER_NO_FULL_COMPLEMENT_SELECTOR_SOURCE_THEOREM",
        )
        self.assertTrue(OUT.exists())
        self.assertTrue(MD.exists())
        self.assertEqual(len(self.theorem["coverage_rows"]), 2)
        self.assertEqual(self.theorem["total_p2451_interval_complement_cell_count"], 99882)
        self.assertEqual(self.theorem["total_p2456_decimal_boundary_replayed_cell_count"], 36)
        self.assertEqual(self.theorem["total_unreplayed_by_decimal_boundary_chain_cell_count"], 99846)

    def test_gap_is_positive_and_not_full_coverage(self) -> None:
        self.assertTrue(self.theorem["p2458_all_weakest_cells_aligned_inherited"])
        self.assertTrue(self.theorem["coverage_gap_positive"])
        self.assertFalse(self.theorem["decimal_boundary_chain_is_full_complement_coverage"])
        self.assertGreater(Decimal(self.theorem["total_unreplayed_coverage_gap_ratio"]), Decimal("0.99"))
        self.assertLess(Decimal(self.theorem["total_decimal_boundary_coverage_ratio"]), Decimal("0.001"))
        for row in self.theorem["coverage_rows"]:
            self.assertFalse(row["decimal_boundary_chain_is_full_family_complement_coverage"])
            self.assertGreater(row["unreplayed_by_decimal_boundary_chain_cell_count"], 0)

    def test_hard_limits(self) -> None:
        self.assertTrue(self.theorem["coverage_gap_ledger_exported"])
        self.assertFalse(self.theorem["directed_rounding_interval_theorem_exported_by_this_certificate"])
        self.assertFalse(self.theorem["symbolic_root_exclusion_theorem_exported_by_this_certificate"])
        self.assertFalse(self.theorem["decimal_full_complement_replay_exported_by_this_certificate"])
        self.assertFalse(self.theorem["pointwise_coordinate_selector_exported_by_this_certificate"])
        self.assertFalse(self.theorem["strict_observable_source_constraint_exported_by_this_certificate"])
        self.assertFalse(self.theorem["gauge_slice_theorem_exported_by_this_certificate"])
        self.assertFalse(self.theorem["strict_physical_value_generator_exported"])
        self.assertFalse(self.theorem["qw2191_discharged"])
        self.assertFalse(self.theorem["role_bearing_ltotal_exported"])
        self.assertFalse(self.theorem["toe_closure_exported"])
        self.assertTrue(all(self.payload["gatekeeper_checks"].values()))

    def test_docs(self) -> None:
        self.assertIn("coverage-gap ledger certificate", MD.read_text(encoding="utf-8"))
        self.assertIn("P2459/S1409", (ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md").read_text(encoding="utf-8"))
        self.assertIn("P2459/S1409", (ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md").read_text(encoding="utf-8"))


if __name__ == "__main__":
    unittest.main()
