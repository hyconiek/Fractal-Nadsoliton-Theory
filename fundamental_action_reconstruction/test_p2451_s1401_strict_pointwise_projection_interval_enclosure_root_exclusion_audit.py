import json
import subprocess
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent
SCRIPT = ROOT / "p2451_s1401_strict_pointwise_projection_interval_enclosure_root_exclusion_audit.py"
OUT = ROOT / "generated" / "p2451_s1401_strict_pointwise_projection_interval_enclosure_root_exclusion_audit.json"
MD = ROOT / "generated" / "p2451_s1401_strict_pointwise_projection_interval_enclosure_root_exclusion_audit.md"
P2450 = ROOT / "generated" / "p2450_s1400_strict_pointwise_projection_root_isolation_margin_certificate.json"


class P2451StrictPointwiseProjectionIntervalEnclosureRootExclusionAuditTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        if not P2450.exists():
            subprocess.run([sys.executable, str(ROOT / "p2450_s1400_strict_pointwise_projection_root_isolation_margin_certificate.py")], check=True)
        subprocess.run([sys.executable, str(SCRIPT)], check=True)
        cls.payload = json.loads(OUT.read_text(encoding="utf-8"))
        cls.cert = cls.payload["strict_pointwise_projection_interval_enclosure_root_exclusion_audit"]
        cls.theorem = cls.cert["theorem_export"]
        cls.zero_projection = cls.theorem["zero_projection_amplitude_interval_audit"]
        cls.stationary_factor = cls.theorem["stationary_factor_interval_audit"]

    def test_identity_outputs_and_cell_counts(self) -> None:
        self.assertEqual(self.payload["packet_id"], "P2451")
        self.assertEqual(self.payload["stage_id"], "S1401")
        self.assertEqual(self.payload["status"], "PASS_STRICT_POINTWISE_PROJECTION_INTERVAL_ENCLOSURE_ROOT_EXCLUSION_AUDIT_NO_EXACT_SELECTOR_SOURCE_THEOREM")
        self.assertTrue(OUT.exists())
        self.assertTrue(MD.exists())
        self.assertGreater(self.zero_projection["cell_count"], 49000)
        self.assertGreater(self.stationary_factor["cell_count"], 49000)

    def test_interval_cells_exclude_zero(self) -> None:
        self.assertEqual(self.zero_projection["zero_containing_cell_count"], 0)
        self.assertEqual(self.stationary_factor["zero_containing_cell_count"], 0)
        self.assertTrue(self.theorem["all_interval_cells_exclude_zero"])
        self.assertGreater(self.zero_projection["minimum_separation_from_zero"], 0.0)
        self.assertGreater(self.stationary_factor["minimum_separation_from_zero"], 0.0)
        self.assertGreater(self.theorem["minimum_interval_separation_from_zero_across_families"], 0.0)

    def test_hard_limits(self) -> None:
        self.assertTrue(self.theorem["floating_interval_enclosure_audit_exported"])
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
        self.assertIn("interval-enclosure root-exclusion audit", MD.read_text(encoding="utf-8"))
        self.assertIn("P2451/S1401", (ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md").read_text(encoding="utf-8"))
        self.assertIn("P2451/S1401", (ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md").read_text(encoding="utf-8"))


if __name__ == "__main__":
    unittest.main()
