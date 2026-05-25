from __future__ import annotations

import json
import subprocess
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent
P1848 = ROOT / "p1848_s798_strict_gravity_componentwise_variation_and_counterterm_witness_checkpoint.py"
P1950 = ROOT / "p1950_s900_strict_renormalization_exact_integration.py"
P2027 = ROOT / "p2027_s977_strict_b1_rank3_gb_null_adaptive_quadrature_witness.py"
SCRIPT = ROOT / "p2028_s978_strict_b1_gb_quotient_counterterm_identifiability_theorem.py"
OUT = ROOT / "generated" / "p2028_s978_strict_b1_gb_quotient_counterterm_identifiability_theorem.json"


class TestP2028S978(unittest.TestCase):
    def test_gb_quotient_counterterm_identifiability_theorem(self) -> None:
        subprocess.run([sys.executable, str(P1848)], check=True)
        subprocess.run([sys.executable, str(P1950)], check=True)
        subprocess.run([sys.executable, str(P2027)], check=True)
        subprocess.run([sys.executable, str(SCRIPT)], check=True)
        data = json.loads(OUT.read_text(encoding="utf-8"))

        self.assertEqual(data["schema_version"], "p2028_s978_v1")
        self.assertEqual(data["status"], "OPEN_OBSTRUCTION_WITH_TRACE")
        self.assertEqual(data["local_verdict"], "PASS_QUOTIENT_THEOREM_WITH_TRACE")
        self.assertEqual(
            data["result_kind"],
            "PASS_B1_SCALAR_GB_QUOTIENT_COUNTERTERM_CLASS_THEOREM__TENSOR_EXTENSION_OPEN",
        )

        exact = data["exact_linear_algebra"]
        self.assertEqual(exact["T_null_vector"], ["0", "0", "0"])
        self.assertTrue(exact["T_section_equals_identity"])
        self.assertTrue(exact["projector_idempotent"])
        self.assertTrue(exact["projector_null_vector_zero"])

        quotient = data["quotient_space"]
        target = quotient["target_quotient_coefficients"]
        self.assertAlmostEqual(target["R2_bar"], 1.0, places=12)
        self.assertAlmostEqual(target["Ric2_bar"], 0.0, places=12)
        self.assertAlmostEqual(target["Riem2_bar"], 0.0, places=12)
        self.assertLessEqual(
            quotient["minimum_norm_quotient_gap_linf"],
            data["tau_family_invariance"]["numerical_tolerance"],
        )

        tau = data["tau_family_invariance"]
        self.assertTrue(tau["tau_family_invariance_pass"])
        self.assertLessEqual(tau["max_quotient_gap_linf"], tau["numerical_tolerance"])
        self.assertLessEqual(tau["max_full_system_residual_linf"], tau["numerical_tolerance"])
        self.assertGreaterEqual(len(tau["tau_samples"]), 6)

        checks = data["gatekeeper_checks"]
        self.assertTrue(checks["p2027_rank3_null1"])
        self.assertTrue(checks["minimum_norm_rep_maps_to_same_quotient"])
        self.assertFalse(checks["independent_a_GB_identified"])
        self.assertFalse(checks["full_tensor_projection_claimed"])
        self.assertFalse(checks["toe_closure_claimed"])
        self.assertIn("quotient-class theorem only", data["false_pass_guard"])
        self.assertIn("tensor-resolved curvature-operator projection", data["theorem_scope"]["not_licensed"])


if __name__ == "__main__":
    unittest.main()
