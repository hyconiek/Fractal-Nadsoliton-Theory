import json
import subprocess
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent
SCRIPT = ROOT / "p2442_s1392_strict_moment_coefficient_local_identifiability_nullspace_certificate.py"
OUT = ROOT / "generated" / "p2442_s1392_strict_moment_coefficient_local_identifiability_nullspace_certificate.json"
MD = ROOT / "generated" / "p2442_s1392_strict_moment_coefficient_local_identifiability_nullspace_certificate.md"
P2441 = ROOT / "generated" / "p2441_s1391_strict_moment_coefficient_phase_sensitivity_rank_certificate.json"


class P2442StrictMomentCoefficientLocalIdentifiabilityNullspaceCertificateTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        if not P2441.exists():
            subprocess.run([sys.executable, str(ROOT / "p2441_s1391_strict_moment_coefficient_phase_sensitivity_rank_certificate.py")], check=True)
        subprocess.run([sys.executable, str(SCRIPT)], check=True)
        cls.payload = json.loads(OUT.read_text(encoding="utf-8"))
        cls.cert = cls.payload["strict_moment_coefficient_local_identifiability_nullspace_certificate"]
        cls.theorem = cls.cert["theorem_export"]

    def test_identity_outputs_and_nullspace_size(self) -> None:
        self.assertEqual(self.payload["packet_id"], "P2442")
        self.assertEqual(self.payload["stage_id"], "S1392")
        self.assertEqual(self.payload["status"], "PASS_STRICT_MOMENT_COEFFICIENT_LOCAL_NULLSPACE_NO_GENERATOR")
        self.assertTrue(OUT.exists())
        self.assertTrue(MD.exists())
        self.assertEqual(self.theorem["input_jacobian_rank"], 3)
        self.assertEqual(self.theorem["nullspace_dimension"], 1)
        self.assertEqual(len(self.theorem["nullspace_certificate"]["normalized_null_vector"]), 4)

    def test_null_direction_and_perturbation_witness(self) -> None:
        self.assertLess(self.theorem["nullspace_certificate"]["max_abs_linear_residual"], 1e-9)
        self.assertEqual(len(self.theorem["null_perturbation_rows"]), 2)
        self.assertLess(self.theorem["max_abs_coefficient_delta_under_null_perturbations"], 1e-4)
        self.assertGreater(self.theorem["max_abs_kernel_sample_delta_under_null_perturbations"], 1e-5)
        self.assertFalse(self.theorem["moment_map_locally_injective_for_four_strict_parameters"])
        self.assertTrue(self.theorem["extra_constraint_or_gauge_redundancy_theorem_required"])

    def test_hard_limits(self) -> None:
        self.assertFalse(self.theorem["strict_kernel_to_coefficient_map_theorem_exported_by_this_certificate"])
        self.assertFalse(self.theorem["strict_physical_value_generator_exported"])
        self.assertFalse(self.theorem["qw2191_discharged"])
        self.assertFalse(self.theorem["role_bearing_ltotal_exported"])
        self.assertFalse(self.theorem["toe_closure_exported"])
        self.assertTrue(all(self.payload["gatekeeper_checks"].values()))

    def test_docs(self) -> None:
        self.assertIn("local-identifiability nullspace certificate", MD.read_text(encoding="utf-8"))
        self.assertIn("P2442/S1392", (ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md").read_text(encoding="utf-8"))
        self.assertIn("P2442/S1392", (ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md").read_text(encoding="utf-8"))


if __name__ == "__main__":
    unittest.main()
