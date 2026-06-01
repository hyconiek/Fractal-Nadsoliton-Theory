import json
import subprocess
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "fundamental_action_reconstruction" / "scratch" / "bridge_strict_completion_legacy_to_strict_finite_bridge_assembly_certificate_probe.py"
REPORT_JSON = ROOT / "fundamental_action_reconstruction" / "scratch" / "bridge_strict_completion_legacy_to_strict_finite_bridge_assembly_certificate_report.json"
REPORT_MD = ROOT / "fundamental_action_reconstruction" / "scratch" / "bridge_strict_completion_legacy_to_strict_finite_bridge_assembly_certificate_report.md"


class LegacyToStrictFiniteBridgeAssemblyCertificateProbeTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        subprocess.run([sys.executable, str(SCRIPT)], check=True, cwd=ROOT)
        cls.payload = json.loads(REPORT_JSON.read_text(encoding="utf-8"))

    def test_result_kind_sources_and_definition(self):
        payload = self.payload
        self.assertEqual(
            payload["result_kind"],
            "SCRATCH_STRICT_COMPLETION_LEGACY_TO_STRICT_FINITE_BRIDGE_ASSEMBLY_CERTIFICATE__COMPARISON_SCOPE_ONLY",
        )
        self.assertIn("finite-apd-bridge-assembly", payload["status"])
        self.assertEqual(set(payload["source_reports"]), {
            "necessity",
            "amplitude_scalar_normalization",
            "finite_diagonal_completion_map",
            "phase_frequency_affine_transport",
            "damping_parameter_identifiability",
            "component_gap_matrix",
            "legacy_bridge_guardrail",
        })
        self.assertIn("t99_bridge_closure_witness_target_spec", payload["source_specs"])
        self.assertIn("assembled bridge witness", payload["grep_disambiguation"]["searched_terms"])
        self.assertIn("assembles", payload["grep_disambiguation"]["finding"])

        definition = payload["assembly_definition"]
        self.assertIn("K_legacy_ont", definition["legacy_kernel"])
        self.assertIn("K_strict_gate", definition["strict_kernel"])
        self.assertEqual(definition["assembled_completion_factor"], "Q_assembly(d)=A(d)*P(d)*D(d)")
        self.assertIn("K_strict_gate(d)=K_legacy_ont(d)*Q_assembly(d)", definition["finite_identity_checked"])
        self.assertIn("not strict dynamical source theorem", definition["scope"])

    def test_rows_summary_and_cross_checks(self):
        payload = self.payload
        rows = payload["assembly_rows"]
        self.assertEqual(len(rows), 12)
        self.assertEqual([row["d"] for row in rows], list(range(12)))
        self.assertTrue(all(abs(row["assembled_Q_minus_diagonal_Q"]) <= 1e-12 for row in rows))
        self.assertTrue(all(abs(row["assembled_Q_minus_necessity_factor_product"]) <= 1e-12 for row in rows))
        self.assertTrue(all(abs(row["reconstruction_residual"]) <= 1e-12 for row in rows))
        self.assertTrue(all(abs(row["log_abs_Q_minus_log_abs_A_plus_log_abs_P_plus_log_abs_D"]) <= 1e-12 for row in rows))
        self.assertTrue(all(row["phase_sign_bit"] == row["diagonal_sign_bit"] for row in rows))
        self.assertTrue(all(row["phase_bit_matches_prior_affine_report"] for row in rows))

        summary = payload["finite_bridge_assembly_summary"]
        self.assertEqual(summary["domain_size"], 12)
        self.assertTrue(summary["finite_kernel_comparison_witness_exported"])
        self.assertTrue(summary["t99_positive_bridge_closure_target_detected"])
        self.assertTrue(summary["comparison_scope_only"])
        self.assertTrue(summary["amplitude_scalar_normalization_inherited"])
        self.assertTrue(summary["phase_affine_transport_inherited"])
        self.assertTrue(summary["damping_beta_eta_identifiability_inherited"])
        self.assertTrue(summary["finite_diagonal_completion_map_inherited"])
        self.assertTrue(summary["assembled_map_reconstructs_strict_kernel"])
        self.assertTrue(summary["assembled_map_matches_finite_diagonal_certificate"])
        self.assertTrue(summary["assembled_map_matches_necessity_apd_product"])
        self.assertTrue(summary["log_abs_decomposition_additive"])
        self.assertTrue(summary["phase_sign_bits_match_diagonal_sign_bits"])
        self.assertTrue(summary["phase_sign_bits_match_prior_affine_report"])
        self.assertTrue(summary["component_gap_still_not_full_bridge"])
        self.assertTrue(summary["guardrail_role_transfer_required_after_full_bridge"])
        self.assertFalse(summary["strict_dynamic_source_exported"])
        self.assertFalse(summary["legacy_role_transfer_exported"])
        self.assertFalse(summary["selector_source_exported"])
        self.assertFalse(summary["full_bridge_theorem_exported"])
        self.assertFalse(summary["toe_closure_claimed"])

        self.assertTrue(all(payload["cross_checks"].values()))

    def test_proof_and_hard_limits(self):
        payload = self.payload
        proof = payload["proof_certificate"]
        self.assertIn("rg was used", proof["grep_step"])
        self.assertIn("Q_assembly(d)=alpha_geo^{-1}", proof["assembly_step"])
        self.assertIn("cancels alpha_geo", proof["identity_step"])
        self.assertIn("Q=diag(K_strict/K_legacy)", proof["diagonal_step"])
        self.assertIn("log|Q|=log(A)+log|P|+log(D)", proof["log_sign_step"])
        self.assertIn("finite T99-style kernel-comparison witness only", proof["scope_step"])

        hard_limits = "\n".join(payload["hard_limits"])
        self.assertIn("No unqualified raw identity", hard_limits)
        self.assertIn("No strict dynamical derivation", hard_limits)
        self.assertIn("No legacy physical-role transfer", hard_limits)
        self.assertIn("No QW-2191 selector discharge", hard_limits)
        self.assertIn("No ToE closure", hard_limits)
        self.assertTrue(REPORT_MD.exists())


if __name__ == "__main__":
    unittest.main()
