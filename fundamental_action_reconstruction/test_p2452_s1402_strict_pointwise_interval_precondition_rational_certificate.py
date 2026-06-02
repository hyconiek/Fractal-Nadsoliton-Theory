import json
import subprocess
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent
SCRIPT = ROOT / "p2452_s1402_strict_pointwise_interval_precondition_rational_certificate.py"
OUT = ROOT / "generated" / "p2452_s1402_strict_pointwise_interval_precondition_rational_certificate.json"
MD = ROOT / "generated" / "p2452_s1402_strict_pointwise_interval_precondition_rational_certificate.md"
P2451 = ROOT / "generated" / "p2451_s1401_strict_pointwise_projection_interval_enclosure_root_exclusion_audit.json"


class P2452StrictPointwiseIntervalPreconditionRationalCertificateTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        if not P2451.exists():
            subprocess.run([sys.executable, str(ROOT / "p2451_s1401_strict_pointwise_projection_interval_enclosure_root_exclusion_audit.py")], check=True)
        subprocess.run([sys.executable, str(SCRIPT)], check=True)
        cls.payload = json.loads(OUT.read_text(encoding="utf-8"))
        cls.cert = cls.payload["strict_pointwise_interval_precondition_rational_certificate"]
        cls.theorem = cls.cert["theorem_export"]

    def test_identity_outputs_and_phase_band(self) -> None:
        self.assertEqual(self.payload["packet_id"], "P2452")
        self.assertEqual(self.payload["stage_id"], "S1402")
        self.assertEqual(self.payload["status"], "PASS_STRICT_POINTWISE_INTERVAL_PRECONDITION_RATIONAL_CERTIFICATE_NO_EXACT_SELECTOR_SOURCE_THEOREM")
        self.assertTrue(OUT.exists())
        self.assertTrue(MD.exists())
        self.assertEqual(self.theorem["exact_phase_interval"]["lo"]["numerator"], 13)
        self.assertEqual(self.theorem["exact_phase_interval"]["lo"]["denominator"], 80)
        self.assertEqual(self.theorem["exact_phase_interval"]["hi"]["numerator"], 873)
        self.assertEqual(self.theorem["exact_phase_interval"]["hi"]["denominator"], 800)
        self.assertTrue(self.theorem["phase_band_inside_open_0_pi_over_2"])
        self.assertTrue(self.theorem["monotone_sin_cos_precondition_certified"])

    def test_domain_and_inheritance(self) -> None:
        self.assertTrue(self.theorem["zero_projection_log_domain_positive"])
        self.assertTrue(self.theorem["stationary_factor_log_domain_positive"])
        self.assertTrue(self.theorem["denominator_positivity_precondition_certified"])
        self.assertEqual(self.theorem["p2451_zero_projection_interval_cell_count"], 49922)
        self.assertEqual(self.theorem["p2451_stationary_factor_interval_cell_count"], 49960)

    def test_hard_limits(self) -> None:
        self.assertTrue(self.theorem["rational_preconditions_exported"])
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
        self.assertIn("interval-precondition rational certificate", MD.read_text(encoding="utf-8"))
        self.assertIn("P2452/S1402", (ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md").read_text(encoding="utf-8"))
        self.assertIn("P2452/S1402", (ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md").read_text(encoding="utf-8"))


if __name__ == "__main__":
    unittest.main()
