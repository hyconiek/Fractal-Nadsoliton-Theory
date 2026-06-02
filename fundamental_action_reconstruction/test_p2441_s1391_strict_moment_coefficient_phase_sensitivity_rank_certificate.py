import json
import subprocess
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent
SCRIPT = ROOT / "p2441_s1391_strict_moment_coefficient_phase_sensitivity_rank_certificate.py"
OUT = ROOT / "generated" / "p2441_s1391_strict_moment_coefficient_phase_sensitivity_rank_certificate.json"
MD = ROOT / "generated" / "p2441_s1391_strict_moment_coefficient_phase_sensitivity_rank_certificate.md"
P2440 = ROOT / "generated" / "p2440_s1390_current_strict_tuple_coefficient_replay_rank_certificate.json"


class P2441StrictMomentCoefficientPhaseSensitivityRankCertificateTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        if not P2440.exists():
            subprocess.run([sys.executable, str(ROOT / "p2440_s1390_current_strict_tuple_coefficient_replay_rank_certificate.py")], check=True)
        subprocess.run([sys.executable, str(SCRIPT)], check=True)
        cls.payload = json.loads(OUT.read_text(encoding="utf-8"))
        cls.cert = cls.payload["strict_moment_coefficient_phase_sensitivity_rank_certificate"]
        cls.theorem = cls.cert["theorem_export"]

    def test_identity_outputs_and_dimensions(self) -> None:
        self.assertEqual(self.payload["packet_id"], "P2441")
        self.assertEqual(self.payload["stage_id"], "S1391")
        self.assertEqual(self.payload["status"], "PASS_STRICT_MOMENT_COEFFICIENT_PHASE_SENSITIVE_NO_GENERATOR")
        self.assertTrue(OUT.exists())
        self.assertTrue(MD.exists())
        self.assertEqual(self.theorem["coefficient_order"], ["lambda_sm_eff", "kappa_gr_eff", "epsilon_mix_eff"])
        self.assertEqual(self.theorem["parameter_order"], ["omega", "phi", "beta", "eta"])
        self.assertEqual(len(self.theorem["jacobian_numeric"]), 3)

    def test_phase_sensitivity_rank_and_sweep(self) -> None:
        self.assertEqual(self.theorem["jacobian_real_rank"], 3)
        self.assertTrue(self.theorem["phi_column_nonzero"])
        self.assertGreater(self.theorem["jacobian_column_norms"]["phi"], 1.0)
        self.assertEqual(self.theorem["phi_sweep_row_count"], 4)
        self.assertTrue(self.theorem["every_phi_sweep_row_changes_a_coefficient"])
        self.assertLess(self.theorem["p1562_comparison"]["max_abs_delta"], 2e-4)

    def test_p2440_obstruction_and_hard_limits(self) -> None:
        self.assertTrue(self.theorem["p2440_phase_null_inherited"])
        self.assertTrue(self.theorem["p2440_not_a_replacement_for_phase_sensitive_moment_route"])
        self.assertFalse(self.theorem["strict_phase_invariance_theorem_exported_by_this_certificate"])
        self.assertFalse(self.theorem["strict_physical_value_generator_exported"])
        self.assertFalse(self.theorem["qw2191_discharged"])
        self.assertFalse(self.theorem["role_bearing_ltotal_exported"])
        self.assertFalse(self.theorem["toe_closure_exported"])
        self.assertTrue(all(self.payload["gatekeeper_checks"].values()))

    def test_docs(self) -> None:
        self.assertIn("phase-sensitivity rank certificate", MD.read_text(encoding="utf-8"))
        self.assertIn("P2441/S1391", (ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md").read_text(encoding="utf-8"))
        self.assertIn("P2441/S1391", (ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md").read_text(encoding="utf-8"))


if __name__ == "__main__":
    unittest.main()
