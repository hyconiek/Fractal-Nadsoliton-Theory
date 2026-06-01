from __future__ import annotations

import json
import subprocess
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "fundamental_action_reconstruction" / "scratch" / "bridge_strict_completion_legacy_to_strict_finite_diagonal_completion_map_certificate_probe.py"
REPORT_JSON = ROOT / "fundamental_action_reconstruction" / "scratch" / "bridge_strict_completion_legacy_to_strict_finite_diagonal_completion_map_certificate_report.json"
REPORT_MD = ROOT / "fundamental_action_reconstruction" / "scratch" / "bridge_strict_completion_legacy_to_strict_finite_diagonal_completion_map_certificate_report.md"


class LegacyToStrictFiniteDiagonalCompletionMapCertificateProbeTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        subprocess.run([sys.executable, str(SCRIPT)], check=True, cwd=ROOT)
        cls.payload = json.loads(REPORT_JSON.read_text(encoding="utf-8"))

    def test_result_kind_sources_and_operator_definition(self):
        payload = self.payload
        self.assertEqual(
            payload["result_kind"],
            "SCRATCH_STRICT_COMPLETION_LEGACY_TO_STRICT_FINITE_DIAGONAL_COMPLETION_MAP_CERTIFICATE__NO_SOURCE_OR_ROLE_TRANSFER",
        )
        self.assertIn("unique-finite-diagonal-completion-map", payload["status"])
        self.assertEqual(set(payload["source_reports"]), {
            "necessity",
            "amplitude_scalar_normalization",
            "damping_compression_separation",
            "positive_factor_sign_separation",
            "component_gap_matrix",
            "legacy_bridge_guardrail",
        })
        self.assertIn("diagonal completion map", payload["grep_disambiguation"]["searched_terms"])
        self.assertIn("unique finite diagonal operator", payload["grep_disambiguation"]["finding"])

        definition = payload["finite_operator_definition"]
        self.assertEqual(definition["domain"], list(range(12)))
        self.assertIn("K_L", definition["legacy_vector_name"])
        self.assertIn("K_S", definition["strict_vector_name"])
        self.assertIn("Q=diag", definition["operator"])
        self.assertEqual(len(definition["diagonal_matrix"]), 12)
        self.assertEqual(len(definition["support_matrix_gf2"]), 12)

    def test_operator_rows_summary_and_cross_checks(self):
        payload = self.payload
        rows = payload["operator_rows"]
        self.assertEqual(len(rows), 12)
        self.assertEqual([row["d"] for row in rows], list(range(12)))
        self.assertTrue(all(abs(row["legacy_value"]) > 1e-14 for row in rows))
        self.assertTrue(all(abs(row["strict_value"]) > 1e-14 for row in rows))
        self.assertTrue(all(abs(row["reconstruction_residual"]) <= 1e-14 for row in rows))
        self.assertTrue(all(abs(row["q_minus_APD"]) <= 1e-14 for row in rows))
        self.assertEqual([row["diagonal_sign_bit"] for row in rows], [0, 0, 1, 1, 1, 1, 0, 0, 1, 1, 0, 0])

        summary = payload["finite_diagonal_completion_summary"]
        self.assertEqual(summary["domain_size"], 12)
        self.assertTrue(summary["legacy_vector_has_no_zero_components"])
        self.assertTrue(summary["strict_vector_has_no_zero_components"])
        self.assertTrue(summary["unique_diagonal_completion_map_exists"])
        self.assertEqual(summary["diagonal_support_rank_over_gf2"], 12)
        self.assertTrue(summary["diagonal_operator_full_rank_on_finite_domain"])
        self.assertNotEqual(summary["diagonal_determinant_nonzero_witness_float"], 0.0)
        self.assertLessEqual(summary["max_abs_reconstruction_residual"], 1e-14)
        self.assertLessEqual(summary["max_abs_q_minus_apd_residual"], 1e-14)
        self.assertTrue(summary["apd_factorization_inherited_exact"])
        self.assertTrue(summary["scalar_only_completion_fails"])
        self.assertGreater(summary["scalar_only_max_abs_residual_using_q0"], 1e-3)
        self.assertTrue(summary["diagonal_sign_bits_match_positive_factor_phase_bits"])
        self.assertTrue(summary["amplitude_scalar_normalization_inherited"])
        self.assertTrue(summary["damping_linear_to_nonlinear_gap_still_open"])
        self.assertTrue(summary["component_gap_matrix_still_not_full_bridge"])
        self.assertFalse(summary["role_transfer_allowed"])
        self.assertFalse(summary["strict_dynamic_source_exported"])
        self.assertFalse(summary["raw_identity_claimed"])
        self.assertFalse(summary["full_bridge_theorem_exported"])

        self.assertTrue(all(payload["cross_checks"].values()))

    def test_proof_and_hard_limits(self):
        payload = self.payload
        proof = payload["proof_certificate"]
        self.assertIn("rg was used", proof["grep_step"])
        self.assertIn("Q=diag(q_d)", proof["existence_step"])
        self.assertIn("unique solution", proof["uniqueness_step"])
        self.assertIn("rank 12", proof["rank_step"])
        self.assertIn("A(d)P(d)D(d)", proof["factorization_step"])
        self.assertIn("not a scalar-only normalization map", proof["non_scalar_step"])
        self.assertIn("not a strict dynamical derivation", proof["theoretical_limit"])

        hard_limits = "\n".join(payload["hard_limits"])
        self.assertIn("No raw identity K_legacy_ont == K_strict_gate", hard_limits)
        self.assertIn("No strict dynamical source", hard_limits)
        self.assertIn("No beta_tors -> chi_11 theorem", hard_limits)
        self.assertIn("No legacy physical-role transfer", hard_limits)
        self.assertIn("No ToE closure", hard_limits)
        self.assertTrue(REPORT_MD.exists())


if __name__ == "__main__":
    unittest.main()
