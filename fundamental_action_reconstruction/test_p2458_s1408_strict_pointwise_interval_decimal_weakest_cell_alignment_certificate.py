import json
import subprocess
import sys
import unittest
from decimal import Decimal
from pathlib import Path

ROOT = Path(__file__).resolve().parent
SCRIPT = ROOT / "p2458_s1408_strict_pointwise_interval_decimal_weakest_cell_alignment_certificate.py"
OUT = ROOT / "generated" / "p2458_s1408_strict_pointwise_interval_decimal_weakest_cell_alignment_certificate.json"
MD = ROOT / "generated" / "p2458_s1408_strict_pointwise_interval_decimal_weakest_cell_alignment_certificate.md"
P2457 = ROOT / "generated" / "p2457_s1407_strict_pointwise_decimal_root_boundary_separation_shape_certificate.json"


class P2458StrictPointwiseIntervalDecimalWeakestCellAlignmentCertificateTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        if not P2457.exists():
            subprocess.run([sys.executable, str(ROOT / "p2457_s1407_strict_pointwise_decimal_root_boundary_separation_shape_certificate.py")], check=True)
        subprocess.run([sys.executable, str(SCRIPT)], check=True)
        cls.payload = json.loads(OUT.read_text(encoding="utf-8"))
        cls.cert = cls.payload["strict_pointwise_interval_decimal_weakest_cell_alignment_certificate"]
        cls.theorem = cls.cert["theorem_export"]

    def test_identity_and_alignment_counts(self) -> None:
        self.assertEqual(self.payload["packet_id"], "P2458")
        self.assertEqual(self.payload["stage_id"], "S1408")
        self.assertEqual(
            self.payload["status"],
            "PASS_STRICT_POINTWISE_INTERVAL_DECIMAL_WEAKEST_CELL_ALIGNMENT_NO_DIRECTED_INTERVAL_SELECTOR_SOURCE_THEOREM",
        )
        self.assertTrue(OUT.exists())
        self.assertTrue(MD.exists())
        self.assertEqual(self.theorem["aligned_family_count"], 2)
        self.assertEqual(len(self.theorem["weakest_cell_alignments"]), 2)

    def test_weakest_cells_match_decimal_boundary_chain(self) -> None:
        self.assertTrue(self.theorem["all_p2451_weakest_cells_match_p2456_boundary_replays"])
        self.assertTrue(self.theorem["all_matched_cells_are_nearest_boundary_cells"])
        self.assertTrue(self.theorem["all_matched_cells_are_covered_by_p2457_shape_audits"])
        self.assertTrue(self.theorem["all_matched_shape_audits_are_strictly_increasing"])
        self.assertTrue(self.theorem["all_matched_shape_audits_are_sign_coherent"])
        self.assertTrue(self.theorem["all_decimal_separations_positive"])
        self.assertGreaterEqual(Decimal(self.theorem["minimum_decimal_minus_float_separation"]), Decimal("0"))
        for alignment in self.theorem["weakest_cell_alignments"]:
            self.assertTrue(alignment["matched_by_p2456_boundary_replay"])
            self.assertEqual(alignment["matched_boundary_row_index"], 0)
            self.assertTrue(alignment["matched_cell_is_nearest_boundary_cell"])
            self.assertTrue(alignment["p2457_shape_audit_found"])

    def test_hard_limits(self) -> None:
        self.assertTrue(self.theorem["backend_chain_alignment_certificate_exported"])
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
        self.assertIn("interval-Decimal weakest-cell alignment certificate", MD.read_text(encoding="utf-8"))
        self.assertIn("P2458/S1408", (ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md").read_text(encoding="utf-8"))
        self.assertIn("P2458/S1408", (ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md").read_text(encoding="utf-8"))


if __name__ == "__main__":
    unittest.main()
