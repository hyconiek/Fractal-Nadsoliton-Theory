from __future__ import annotations

import json
import subprocess
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "fundamental_action_reconstruction" / "scratch" / "bridge_strict_kernel_completion_necessity_certificate_probe.py"
REPORT_JSON = ROOT / "fundamental_action_reconstruction" / "scratch" / "bridge_strict_kernel_completion_necessity_certificate_report.json"
REPORT_MD = ROOT / "fundamental_action_reconstruction" / "scratch" / "bridge_strict_kernel_completion_necessity_certificate_report.md"


class StrictKernelCompletionNecessityCertificateProbeTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        subprocess.run([sys.executable, str(SCRIPT)], check=True, cwd=ROOT)
        cls.payload = json.loads(REPORT_JSON.read_text(encoding="utf-8"))

    def test_result_kind_constants_and_factor_definitions(self):
        payload = self.payload
        self.assertEqual(
            payload["result_kind"],
            "SCRATCH_STRICT_KERNEL_COMPLETION_NECESSITY_CERTIFICATE__FINITE_Z12_NO_DYNAMICAL_DERIVATION",
        )
        self.assertIn("all-three-explicit-factors", payload["status"])
        constants = payload["constants"]
        self.assertAlmostEqual(constants["alpha_geo"], 2.772588722239781)
        self.assertEqual(constants["domain_d_values"], list(range(12)))
        self.assertEqual(constants["strict"]["eta"], 9.0 / 5.0)
        definitions = payload["completion_factor_definitions"]
        self.assertIn("1/alpha_geo", definitions["alpha_normalization"])
        self.assertIn("cos(omega_S", definitions["phase_frequency_transport"])
        self.assertIn("beta_tors*d", definitions["damping_compression"])
        self.assertIn("K_strict_gate", definitions["identity"])

    def test_subset_exhaustion_and_scalar_repair(self):
        summary = self.payload["necessity_summary"]
        self.assertEqual(
            summary["exact_subsets_without_extra_scalar"],
            ["alpha_normalization+phase_frequency_transport+damping_compression"],
        )
        self.assertEqual(
            summary["exact_subsets_up_to_one_global_scalar"],
            [
                "phase_frequency_transport+damping_compression",
                "alpha_normalization+phase_frequency_transport+damping_compression",
            ],
        )
        self.assertEqual(summary["phase_missing_sign_mismatch_union"], [2, 3, 4, 5, 8, 9])
        self.assertGreater(summary["minimum_best_scalar_l2_residual_when_phase_missing"], 0.5)
        self.assertGreater(summary["minimum_best_scalar_l2_residual_when_damping_missing"], 0.5)
        self.assertLess(summary["max_abs_quotient_identity_residual"], 1e-15)
        self.assertTrue(summary["damping_factor_positive"])
        self.assertTrue(summary["damping_factor_strictly_decreasing_after_d0"])
        self.assertIn("global scalar", summary["alpha_is_global_scalar_not_shape_factor"])
        self.assertIn("signs disagree", summary["phase_is_sign_shape_critical"])
        self.assertIn("tail envelope", summary["damping_is_envelope_shape_critical"])

    def test_subset_rows_and_pointwise_quotient_certificate(self):
        subset_rows = self.payload["subset_enumeration"]
        self.assertEqual(len(subset_rows), 8)
        by_label = {row["subset_label"]: row for row in subset_rows}
        full = by_label["alpha_normalization+phase_frequency_transport+damping_compression"]
        self.assertTrue(full["exact_without_extra_scalar"])
        self.assertLess(full["max_abs_residual"], 1e-15)
        phase_damping = by_label["phase_frequency_transport+damping_compression"]
        self.assertFalse(phase_damping["exact_without_extra_scalar"])
        self.assertTrue(phase_damping["exact_up_to_one_global_scalar"])
        self.assertAlmostEqual(phase_damping["best_scalar_repair"]["best_scalar"], 0.36067376022224085)
        no_phase = by_label["alpha_normalization+damping_compression"]
        self.assertEqual(no_phase["sign_mismatch_d_values"], [2, 3, 4, 5, 8, 9])
        self.assertFalse(no_phase["exact_up_to_one_global_scalar"])

        point_rows = self.payload["pointwise_quotient_certificate"]
        self.assertEqual(len(point_rows), 12)
        self.assertEqual(point_rows[0]["d"], 0)
        self.assertEqual(point_rows[-1]["d"], 11)
        self.assertAlmostEqual(point_rows[0]["damping_compression_factor"], 1.0)
        self.assertAlmostEqual(point_rows[-1]["damping_compression_factor"], 0.014623674671644927)
        self.assertGreater(point_rows[-1]["tail_amplification_if_damping_omitted"], 68.0)
        for row in point_rows:
            self.assertAlmostEqual(row["strict_over_legacy_quotient"], row["factor_product"], delta=1e-15)

    def test_blocker_context_proof_guardrails_and_markdown(self):
        blockers = self.payload["blocker_context"]
        self.assertIn("legacy-carrier-completed", blockers["ladder_status"])
        self.assertIn("minimal-premise-lattice", blockers["blocker_lattice_status"])
        self.assertIn("inverse", blockers["inverse_branch_status"])
        self.assertIn("strict_transport_derivation", blockers["what_remains_open"])
        self.assertIn("orientation_chi11_source", blockers["what_remains_open"])
        self.assertIn("role_transfer_theorem", blockers["what_remains_open"])

        proof = self.payload["proof_certificate"]
        self.assertIn("K_strict_gate", proof["quotient_step"])
        self.assertIn("All 2^3 subsets", proof["subset_exhaustion"])
        self.assertIn("missing-alpha", proof["scalar_repair_step"])
        self.assertIn("six legacy-vs-strict sign mismatches", proof["sign_step"])
        self.assertIn("not another factorization-only", proof["nonduplication"])

        hard_limits = "\n".join(self.payload["hard_limits"])
        self.assertIn("K_strict_gate remains the current live/full", hard_limits)
        self.assertIn("No unqualified identity K_legacy_ont == K_strict_gate", hard_limits)
        self.assertIn("No proof derives A(d), P(d), or D(d)", hard_limits)
        self.assertIn("No beta_tors -> chi_11 theorem", hard_limits)
        self.assertIn("No QW-2191 selector discharge", hard_limits)
        self.assertIn("No ToE closure", hard_limits)
        self.assertTrue(REPORT_MD.exists())


if __name__ == "__main__":
    unittest.main()
