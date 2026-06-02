import json
import subprocess
import sys
import unittest
from decimal import Decimal
from pathlib import Path

ROOT = Path(__file__).resolve().parent
SCRIPT = ROOT / "p2453_s1403_strict_pointwise_directed_decimal_weakest_cell_replay_certificate.py"
OUT = ROOT / "generated" / "p2453_s1403_strict_pointwise_directed_decimal_weakest_cell_replay_certificate.json"
MD = ROOT / "generated" / "p2453_s1403_strict_pointwise_directed_decimal_weakest_cell_replay_certificate.md"
P2452 = ROOT / "generated" / "p2452_s1402_strict_pointwise_interval_precondition_rational_certificate.json"


class P2453StrictPointwiseDirectedDecimalWeakestCellReplayCertificateTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        if not P2452.exists():
            subprocess.run([sys.executable, str(ROOT / "p2452_s1402_strict_pointwise_interval_precondition_rational_certificate.py")], check=True)
        subprocess.run([sys.executable, str(SCRIPT)], check=True)
        cls.payload = json.loads(OUT.read_text(encoding="utf-8"))
        cls.cert = cls.payload["strict_pointwise_directed_decimal_weakest_cell_replay_certificate"]
        cls.theorem = cls.cert["theorem_export"]

    def test_identity_and_backend_parameters(self) -> None:
        self.assertEqual(self.payload["packet_id"], "P2453")
        self.assertEqual(self.payload["stage_id"], "S1403")
        self.assertEqual(self.payload["status"], "PASS_STRICT_POINTWISE_DIRECTED_DECIMAL_WEAKEST_CELL_REPLAY_NO_FULL_INTERVAL_SELECTOR_SOURCE_THEOREM")
        self.assertTrue(OUT.exists())
        self.assertTrue(MD.exists())
        self.assertEqual(self.theorem["decimal_precision"], 90)
        self.assertTrue(self.theorem["p2452_phase_precondition_inherited"])

    def test_weakest_cells_exclude_zero(self) -> None:
        zero = self.theorem["zero_projection_weakest_cell_replay"]
        stationary = self.theorem["stationary_factor_weakest_cell_replay"]
        self.assertTrue(zero["decimal_interval_excludes_zero"])
        self.assertTrue(stationary["decimal_interval_excludes_zero"])
        self.assertGreater(Decimal(zero["decimal_separation_from_zero"]), Decimal("0"))
        self.assertGreater(Decimal(stationary["decimal_separation_from_zero"]), Decimal("0"))
        self.assertTrue(self.theorem["both_weakest_cells_exclude_zero_under_decimal_taylor_replay"])
        self.assertTrue(self.theorem["both_weakest_cells_keep_p2451_positive_sign"])

    def test_hard_limits(self) -> None:
        self.assertFalse(self.theorem["full_complement_directed_rounding_interval_theorem_exported_by_this_certificate"])
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
        self.assertIn("directed-decimal weakest-cell replay", MD.read_text(encoding="utf-8"))
        self.assertIn("P2453/S1403", (ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md").read_text(encoding="utf-8"))
        self.assertIn("P2453/S1403", (ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md").read_text(encoding="utf-8"))


if __name__ == "__main__":
    unittest.main()
