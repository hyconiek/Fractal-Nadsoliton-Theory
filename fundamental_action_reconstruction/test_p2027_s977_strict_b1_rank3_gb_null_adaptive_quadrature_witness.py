from __future__ import annotations

import json
import subprocess
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent
P1848 = ROOT / "p1848_s798_strict_gravity_componentwise_variation_and_counterterm_witness_checkpoint.py"
SCRIPT = ROOT / "p2027_s977_strict_b1_rank3_gb_null_adaptive_quadrature_witness.py"
OUT = ROOT / "generated" / "p2027_s977_strict_b1_rank3_gb_null_adaptive_quadrature_witness.json"


class TestP2027S977(unittest.TestCase):
    def test_rank3_gb_null_adaptive_quadrature_witness(self) -> None:
        subprocess.run([sys.executable, str(P1848)], check=True)
        subprocess.run([sys.executable, str(SCRIPT)], check=True)
        data = json.loads(OUT.read_text(encoding="utf-8"))

        self.assertEqual(data["schema_version"], "p2027_s977_v1")
        self.assertEqual(data["status"], "OPEN_OBSTRUCTION_WITH_TRACE")
        self.assertEqual(data["local_verdict"], "PASS_RANK3_QUOTIENT_IDENTIFIABLE_ON_SCALAR_B1_WITH_GB_NULL_TRACE")
        self.assertEqual(data["basis_decomposition"]["full_gram_rank"], 3)
        self.assertEqual(data["basis_decomposition"]["full_gram_nullity"], 1)
        self.assertTrue(data["basis_decomposition"]["symbolic_gb_relation_zero"])
        self.assertGreater(data["basis_decomposition"]["null_vector_cosine_alignment_with_exact"], 0.999999999)

        quotient = data["quotient_identifiability"]
        coeffs = quotient["rank3_quotient_coefficients"]
        self.assertTrue(quotient["quotient_gate_pass"])
        self.assertAlmostEqual(coeffs["a_R2"], 1.0, places=12)
        self.assertAlmostEqual(coeffs["a_Ric2"], 0.0, places=12)
        self.assertAlmostEqual(coeffs["a_Riem2"], 0.0, places=12)
        self.assertAlmostEqual(coeffs["a_GB_representative"], 0.0, places=12)

        quad = data["adaptive_quadrature"]
        self.assertTrue(quad["adaptive_transform_gate_pass"])
        self.assertLessEqual(
            quad["adaptive_relative_gap_primary_vs_check"],
            quad["adaptive_relative_stability_tolerance"],
        )
        transforms = {row["transform"] for row in quad["independent_matrix_entries_directly_integrated"]}
        self.assertIn("left_endpoint_power_5", transforms)

        self.assertFalse(data["full_four_channel_min_norm_family"]["gb_coefficient_identified"])
        self.assertTrue(data["gatekeeper_checks"]["no_four_channel_closure_claimed"])
        self.assertIn("not full renormalization closure", data["false_pass_guard"])
        self.assertIn("full tensor projection closure", data["theorem_scope"]["not_licensed"])


if __name__ == "__main__":
    unittest.main()
