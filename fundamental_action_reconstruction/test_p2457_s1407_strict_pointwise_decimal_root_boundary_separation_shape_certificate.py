import json
import subprocess
import sys
import unittest
from decimal import Decimal
from pathlib import Path

ROOT = Path(__file__).resolve().parent
SCRIPT = ROOT / "p2457_s1407_strict_pointwise_decimal_root_boundary_separation_shape_certificate.py"
OUT = ROOT / "generated" / "p2457_s1407_strict_pointwise_decimal_root_boundary_separation_shape_certificate.json"
MD = ROOT / "generated" / "p2457_s1407_strict_pointwise_decimal_root_boundary_separation_shape_certificate.md"
P2456 = ROOT / "generated" / "p2456_s1406_strict_pointwise_decimal_root_window_boundary_band_replay_certificate.json"


class P2457StrictPointwiseDecimalRootBoundarySeparationShapeCertificateTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        if not P2456.exists():
            subprocess.run([sys.executable, str(ROOT / "p2456_s1406_strict_pointwise_decimal_root_window_boundary_band_replay_certificate.py")], check=True)
        subprocess.run([sys.executable, str(SCRIPT)], check=True)
        cls.payload = json.loads(OUT.read_text(encoding="utf-8"))
        cls.cert = cls.payload["strict_pointwise_decimal_root_boundary_separation_shape_certificate"]
        cls.theorem = cls.cert["theorem_export"]

    def test_identity_and_shape_counts(self) -> None:
        self.assertEqual(self.payload["packet_id"], "P2457")
        self.assertEqual(self.payload["stage_id"], "S1407")
        self.assertEqual(
            self.payload["status"],
            "PASS_STRICT_POINTWISE_DECIMAL_ROOT_BOUNDARY_SEPARATION_SHAPE_NO_FULL_INTERVAL_SELECTOR_SOURCE_THEOREM",
        )
        self.assertTrue(OUT.exists())
        self.assertTrue(MD.exists())
        self.assertTrue(self.theorem["p2456_all_boundary_bands_exclude_zero_inherited"])
        self.assertEqual(self.theorem["total_shape_audited_boundary_band_count"], 6)
        self.assertEqual(self.theorem["total_shape_audited_boundary_cell_count"], 36)

    def test_shape_and_sign_witnesses(self) -> None:
        self.assertEqual(len(self.theorem["boundary_band_shape_audits"]), 6)
        self.assertEqual(len(self.theorem["two_sided_root_window_sign_pairs"]), 3)
        self.assertTrue(self.theorem["all_boundary_bands_strictly_increase_separation_away_from_root_window"])
        self.assertTrue(self.theorem["all_boundary_bands_are_sign_coherent"])
        self.assertTrue(self.theorem["all_two_sided_root_windows_have_opposite_boundary_signs"])
        self.assertGreater(Decimal(self.theorem["minimum_first_separation_difference_across_boundary_bands"]), Decimal("0"))
        for audit in self.theorem["boundary_band_shape_audits"]:
            self.assertTrue(audit["strictly_increasing_separation_away_from_root_window"])
            self.assertTrue(audit["sign_coherent_zero_excluding_band"])
            self.assertEqual(len(audit["separations_ordered_near_to_far"]), 6)
            self.assertEqual(len(audit["first_separation_differences_near_to_far"]), 5)

    def test_hard_limits(self) -> None:
        self.assertTrue(self.theorem["bounded_root_boundary_shape_audit_exported"])
        self.assertFalse(self.theorem["full_complement_directed_rounding_interval_theorem_exported_by_this_certificate"])
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
        self.assertIn("root-boundary separation-shape certificate", MD.read_text(encoding="utf-8"))
        self.assertIn("P2457/S1407", (ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md").read_text(encoding="utf-8"))
        self.assertIn("P2457/S1407", (ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md").read_text(encoding="utf-8"))


if __name__ == "__main__":
    unittest.main()
