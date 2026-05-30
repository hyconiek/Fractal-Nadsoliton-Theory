from __future__ import annotations

import json
import subprocess
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "fundamental_action_reconstruction" / "scratch" / "bridge_strict_completion_damping_continuous_monotonicity_certificate_probe.py"
REPORT_JSON = ROOT / "fundamental_action_reconstruction" / "scratch" / "bridge_strict_completion_damping_continuous_monotonicity_certificate_report.json"
REPORT_MD = ROOT / "fundamental_action_reconstruction" / "scratch" / "bridge_strict_completion_damping_continuous_monotonicity_certificate_report.md"


class StrictCompletionDampingContinuousMonotonicityCertificateProbeTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        subprocess.run([sys.executable, str(SCRIPT)], check=True, cwd=ROOT)
        cls.payload = json.loads(REPORT_JSON.read_text(encoding="utf-8"))

    def test_result_kind_constants_and_derivative_certificate(self):
        payload = self.payload
        self.assertEqual(
            payload["result_kind"],
            "SCRATCH_STRICT_COMPLETION_DAMPING_CONTINUOUS_MONOTONICITY_CERTIFICATE__CALCULUS_ENVELOPE_ONLY",
        )
        self.assertIn("positive-and-strictly-decreasing", payload["status"])
        constants = payload["constants"]
        self.assertAlmostEqual(constants["beta_tors"], 0.01)
        self.assertAlmostEqual(constants["beta"], 1.0)
        self.assertAlmostEqual(constants["eta"], 9.0 / 5.0)
        self.assertEqual(constants["certified_continuous_interval"], [1.0, 11.0])
        cert = payload["derivative_certificate"]
        self.assertIn("D'(x)=N(x)", cert["derivative_formula"])
        self.assertIn("beta_tors - eta", cert["inequality_for_x_ge_1"])
        self.assertAlmostEqual(cert["numerator_upper_bound_for_x_ge_1"], -1.79)
        self.assertIn("D'(x)<0", cert["conclusion"])

    def test_monotonicity_summary_grid_and_edges(self):
        summary = self.payload["monotonicity_summary"]
        self.assertAlmostEqual(summary["D_1"], 0.505)
        self.assertAlmostEqual(summary["D_11"], 0.014623674671644927)
        self.assertAlmostEqual(summary["D_1_over_D_11"], 34.533043939987735)
        self.assertTrue(summary["continuous_positive_certificate"])
        self.assertTrue(summary["continuous_strictly_decreasing_certificate"])
        self.assertTrue(summary["all_grid_derivative_numerators_negative"])
        self.assertTrue(summary["all_grid_derivatives_negative"])
        self.assertTrue(summary["all_integer_edge_drops_positive"])
        self.assertLess(summary["max_grid_derivative_numerator"], 0.0)
        self.assertGreater(summary["min_integer_edge_drop"], 0.0)

        grid = self.payload["grid_derivative_rows"]
        self.assertEqual(len(grid), 11)
        for row in grid:
            self.assertTrue(row["derivative_denominator_positive"])
            self.assertTrue(row["N_x_negative"])
            self.assertTrue(row["D_prime_x_negative"])

        edges = self.payload["integer_edge_drop_rows"]
        self.assertEqual(len(edges), 10)
        for row in edges:
            self.assertGreater(row["drop"], 0.0)
            self.assertGreater(row["ratio_left_over_right"], 1.0)

    def test_blocker_context_proof_guardrails_and_markdown(self):
        blockers = self.payload["blocker_context"]
        self.assertIn("all-three-explicit-factors", blockers["necessity_status"])
        self.assertIn("unique-finite-node-transport", blockers["transport_cocycle_status"])
        self.assertIn("phase-transport-sign-flips", blockers["phase_zero_status"])
        self.assertIn("low-order-autonomous", blockers["low_order_status"])
        self.assertIn("still_open", blockers["strict_damping_derivation_status"])
        self.assertIn("orientation_chi11_source", blockers["still_open"])
        self.assertIn("role_transfer_theorem", blockers["still_open"])

        proof = self.payload["proof_certificate"]
        self.assertIn("D(x)>0", proof["positivity_step"])
        self.assertIn("Differentiating D(x)", proof["derivative_step"])
        self.assertIn("N(x)<=beta_tors-eta", proof["negative_numerator_step"])
        self.assertIn("cannot create the phase sign flips", proof["phase_separation_step"])
        self.assertIn("not another factor-subset", proof["nonduplication"])

        hard_limits = "\n".join(self.payload["hard_limits"])
        self.assertIn("K_strict_gate remains the current live/full", hard_limits)
        self.assertIn("No unqualified identity K_legacy_ont == K_strict_gate", hard_limits)
        self.assertIn("No proof derives beta, eta, or D(x)", hard_limits)
        self.assertIn("No beta_tors -> chi_11 theorem", hard_limits)
        self.assertIn("No QW-2191 selector discharge", hard_limits)
        self.assertIn("No ToE closure", hard_limits)
        self.assertTrue(REPORT_MD.exists())


if __name__ == "__main__":
    unittest.main()
