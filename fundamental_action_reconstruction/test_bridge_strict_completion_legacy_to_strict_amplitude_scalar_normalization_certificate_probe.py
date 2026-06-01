from __future__ import annotations

import json
import subprocess
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "fundamental_action_reconstruction" / "scratch" / "bridge_strict_completion_legacy_to_strict_amplitude_scalar_normalization_certificate_probe.py"
REPORT_JSON = ROOT / "fundamental_action_reconstruction" / "scratch" / "bridge_strict_completion_legacy_to_strict_amplitude_scalar_normalization_certificate_report.json"
REPORT_MD = ROOT / "fundamental_action_reconstruction" / "scratch" / "bridge_strict_completion_legacy_to_strict_amplitude_scalar_normalization_certificate_report.md"


class LegacyToStrictAmplitudeScalarNormalizationCertificateProbeTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        subprocess.run([sys.executable, str(SCRIPT)], check=True, cwd=ROOT)
        cls.payload = json.loads(REPORT_JSON.read_text(encoding="utf-8"))

    def test_result_kind_sources_and_definition(self):
        payload = self.payload
        self.assertEqual(
            payload["result_kind"],
            "SCRATCH_STRICT_COMPLETION_LEGACY_TO_STRICT_AMPLITUDE_SCALAR_NORMALIZATION_CERTIFICATE__NO_ROLE_TRANSFER",
        )
        self.assertIn("legacy-alpha-geo-factors-out", payload["status"])
        self.assertEqual(set(payload["source_reports"]), {
            "agents_guardrail",
            "s2_priority_packet",
            "alpha_geo_strict_derived",
            "necessity",
            "component_gap_matrix",
            "legacy_bridge_guardrail",
        })
        searched = "\n".join(payload["grep_disambiguation"]["searched_terms"])
        self.assertIn("amplitude normalization", searched)
        self.assertIn("amplitude nonabsorption", searched)
        self.assertIn("narrower", payload["grep_disambiguation"]["finding"])

        definition = payload["normalization_definition"]
        self.assertIn("alpha_geo", definition["legacy_kernel"])
        self.assertIn("alpha_geo^{-1}", definition["legacy_shape_after_scalar_normalization"])
        self.assertEqual(definition["alpha_geo"], "4 ln(2)")
        self.assertEqual(definition["domain"], list(range(12)))
        self.assertIn("3d == 4 mod 12", definition["cos_zero_congruence"])

    def test_rows_summary_and_cross_checks(self):
        payload = self.payload
        rows = payload["node_rows"]
        self.assertEqual(len(rows), 12)
        self.assertEqual([row["d"] for row in rows], list(range(12)))
        self.assertTrue(all(row["legacy_denominator"] > 0 for row in rows))
        self.assertTrue(all(row["legacy_shape_L_d"] != 0 for row in rows))
        self.assertTrue(all(abs(row["normalization_residual"]) < 1e-15 for row in rows))
        self.assertTrue(all(row["sign_preserved_by_positive_alpha"] for row in rows))

        summary = payload["amplitude_scalar_normalization_summary"]
        self.assertTrue(summary["legacy_alpha_geo_visible"])
        self.assertTrue(summary["strict_alpha_geo_source_loaded"])
        self.assertTrue(summary["legacy_shape_nonzero_on_domain"])
        self.assertTrue(summary["legacy_denominator_positive_on_domain"])
        self.assertTrue(summary["cos_zero_congruence_has_no_domain_solution"])
        self.assertTrue(summary["alpha_geo_positive"])
        self.assertTrue(summary["alpha_inverse_normalization_residual_zero_formally"])
        self.assertLess(summary["alpha_inverse_normalization_residual_max_abs_float"], 1e-15)
        self.assertTrue(summary["ratio_K_over_shape_constant_alpha_geo_formally"])
        self.assertLess(summary["ratio_K_over_shape_max_float_deviation"], 1e-15)
        self.assertTrue(summary["sign_pattern_preserved_by_positive_alpha"])
        self.assertTrue(summary["necessity_exact_apd_inherited"])
        self.assertTrue(summary["component_gap_amplitude_row_present"])
        self.assertTrue(summary["scalar_normalization_witness_exported"])
        self.assertFalse(summary["full_strict_A_d_derivation_exported"])
        self.assertFalse(summary["legacy_role_transfer_allowed"])
        self.assertFalse(summary["raw_kernel_identity_claimed"])
        self.assertFalse(summary["full_bridge_theorem_exported"])

        self.assertTrue(all(payload["cross_checks"].values()))

    def test_proof_and_hard_limits(self):
        payload = self.payload
        proof = payload["proof_certificate"]
        self.assertIn("rg was used", proof["grep_step"])
        self.assertIn("K_legacy_ont(d)=alpha_geo*L_shape(d)", proof["factorization_step"])
        self.assertIn("3d==4 mod 12 has no solution", proof["nonzero_step"])
        self.assertIn("preserves the legacy sign pattern", proof["positivity_step"])
        self.assertIn("alpha_geo_strict_derived_v1=4 ln(2)", proof["strict_source_step"])
        self.assertIn("amplitude row now has an explicit scalar-normalization witness", proof["bridge_meaning_step"])
        self.assertIn("not a transfer", proof["theoretical_limit"])

        hard_limits = "\n".join(payload["hard_limits"])
        self.assertIn("No raw identity K_legacy_ont == K_strict_gate", hard_limits)
        self.assertIn("No full strict A(d) dynamical derivation", hard_limits)
        self.assertIn("No legacy alpha_geo physical-role transfer", hard_limits)
        self.assertIn("No QW-2191 selector discharge", hard_limits)
        self.assertTrue(REPORT_MD.exists())


if __name__ == "__main__":
    unittest.main()
