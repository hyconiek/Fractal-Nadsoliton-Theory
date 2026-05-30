from __future__ import annotations

import json
import subprocess
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "fundamental_action_reconstruction" / "scratch" / "bridge_completed_strict_kernel_factorization_certificate_probe.py"
REPORT_JSON = ROOT / "fundamental_action_reconstruction" / "scratch" / "bridge_completed_strict_kernel_factorization_certificate_report.json"
REPORT_MD = ROOT / "fundamental_action_reconstruction" / "scratch" / "bridge_completed_strict_kernel_factorization_certificate_report.md"


class CompletedStrictKernelFactorizationCertificateProbeTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        subprocess.run([sys.executable, str(SCRIPT)], check=True, cwd=ROOT)
        cls.payload = json.loads(REPORT_JSON.read_text(encoding="utf-8"))

    def test_result_kind_kernel_reading_and_constants(self):
        payload = self.payload
        self.assertEqual(
            payload["result_kind"],
            "SCRATCH_COMPLETED_STRICT_KERNEL_FACTORIZATION_CERTIFICATE__EXACT_FACTOR_TARGET_NOT_BRIDGE_THEOREM",
        )
        self.assertEqual(payload["status"], "strict-kernel-current-full-form-with-explicit-legacy-completion-factorization")
        reading = payload["kernel_reading"]
        self.assertIn("K_strict_gate", reading["current_live_kernel"])
        self.assertIn("historical/nadsoliton-characteristic carrier", reading["legacy_reading"])
        self.assertIn("explicit completion factors", reading["completion_reading"])
        self.assertIn("Pointwise factorization is exact", reading["guarded_identity"])
        constants = payload["constants"]
        self.assertAlmostEqual(constants["alpha_geo"], 2.772588722239781)
        self.assertEqual(constants["target_q_power"], "256/243")
        self.assertEqual(constants["target_eta"], "9/5")
        self.assertEqual(constants["balanced_ledger"], [2, 2, 2, 1, 1])

    def test_factorization_summary_and_rows(self):
        summary = self.payload["factorization_summary"]
        self.assertEqual(summary["domain_d_values"], list(range(1, 12)))
        self.assertLess(summary["max_abs_residual_from_legacy_full"], 1e-15)
        self.assertLess(summary["max_abs_residual_from_legacy_reduced"], 1e-15)
        self.assertTrue(summary["residual_tolerance_pass"])
        self.assertTrue(summary["damping_compression_factor_positive"])
        self.assertTrue(summary["damping_compression_factor_strictly_decreasing"])
        self.assertAlmostEqual(summary["damping_factor_d1"], 0.505)
        self.assertAlmostEqual(summary["damping_factor_d11"], 0.014623674671644927)
        self.assertAlmostEqual(summary["damping_factor_d1_over_d11"], 34.533043939987735)
        self.assertEqual(summary["phase_transport_sign_pattern"], [1, -1, -1, -1, -1, 1, 1, -1, -1, 1, 1])

        rows = self.payload["factorization_rows"]
        self.assertEqual(len(rows), 11)
        for row in rows:
            self.assertAlmostEqual(row["reconstructed_from_legacy_full"], row["strict_kernel"], delta=1e-15)
            self.assertAlmostEqual(row["reconstructed_from_legacy_reduced"], row["strict_kernel"], delta=1e-15)
            self.assertGreater(row["damping_compression_factor"], 0.0)
        self.assertAlmostEqual(rows[0]["amplitude_factor_alpha_removal"], 0.36067376022224085)
        self.assertAlmostEqual(rows[0]["damping_compression_factor"], 0.505)

    def test_upstream_context_and_interpretation(self):
        upstream = self.payload["upstream_bridge_context"]
        self.assertIn("legacy-kernel-restored", upstream["reactivation_status"])
        self.assertEqual(upstream["one_bit_frontier"]["frontier_name"], "T_beta_tors_orientation_exports_chi11_or_nonbridge")
        self.assertIn("no full-Aut invariant source exports chi_11", upstream["reynolds_full_aut_obstruction"])
        self.assertTrue(upstream["reynolds_annihilator_zero"])

        interpretation = self.payload["completion_factor_interpretation"]
        self.assertIn("1/alpha_geo", interpretation["amplitude"])
        self.assertIn("strict gate phase", interpretation["phase_frequency"])
        self.assertIn("upgrading legacy hyperbolic damping", interpretation["damping_compression"])
        self.assertIn("does not supply beta_tors -> chi_11", interpretation["orientation_bit"])

    def test_proof_guardrails_and_markdown(self):
        proof = self.payload["exact_proof_certificate"]
        self.assertIn("K_strict(d) - C_full(d)*K_legacy_full(d) = 0", proof["factor_identity"])
        self.assertIn("explicit completion factors", proof["completion_objects"])
        self.assertIn("strict completion of the legacy carrier", proof["computational_positive"])
        self.assertIn("not the beta_tors->chi_11 theorem", proof["theoretical_limit"])

        ontology = self.payload["ontology_guardrail"]
        self.assertIn("nadsoliton itself", ontology["allowed_reading"])
        self.assertIn("No separate informational layer", ontology["forbidden_reading"])

        hard_limits = "\n".join(self.payload["hard_limits"])
        self.assertIn("K_strict_gate is the current full form", hard_limits)
        self.assertIn("No unqualified identity K_legacy_ont == K_strict_gate", hard_limits)
        self.assertIn("No beta_tors -> chi_11 theorem", hard_limits)
        self.assertIn("No proof derives the completion factors", hard_limits)
        self.assertIn("No legacy physical-role transfer", hard_limits)
        self.assertIn("No QW-2191 discharge", hard_limits)
        self.assertIn("No ToE closure", hard_limits)
        self.assertTrue(REPORT_MD.exists())


if __name__ == "__main__":
    unittest.main()
