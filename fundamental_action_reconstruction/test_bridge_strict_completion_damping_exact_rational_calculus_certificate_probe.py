from __future__ import annotations

import json
import subprocess
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "fundamental_action_reconstruction" / "scratch" / "bridge_strict_completion_damping_exact_rational_calculus_certificate_probe.py"
REPORT_JSON = ROOT / "fundamental_action_reconstruction" / "scratch" / "bridge_strict_completion_damping_exact_rational_calculus_certificate_report.json"
REPORT_MD = ROOT / "fundamental_action_reconstruction" / "scratch" / "bridge_strict_completion_damping_exact_rational_calculus_certificate_report.md"


class StrictCompletionDampingExactRationalCalculusCertificateProbeTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        subprocess.run([sys.executable, str(SCRIPT)], check=True, cwd=ROOT)
        cls.payload = json.loads(REPORT_JSON.read_text(encoding="utf-8"))

    def test_result_kind_sources_and_parameters(self):
        payload = self.payload
        self.assertEqual(
            payload["result_kind"],
            "SCRATCH_STRICT_COMPLETION_DAMPING_EXACT_RATIONAL_CALCULUS_CERTIFICATE__NO_FLOAT_DERIVATIVE_SIGN",
        )
        self.assertIn("exact-rational-bound", payload["status"])
        self.assertEqual(set(payload["source_reports"]), {
            "continuous_damping_monotonicity_certificate",
            "cell_sign_certificate",
            "phase_sign_z2_coboundary_certificate",
        })
        params = payload["rational_parameters"]
        self.assertEqual(params["beta"]["text"], "1/1")
        self.assertEqual(params["eta"]["text"], "9/5")
        self.assertEqual(params["beta_tors"]["text"], "1/100")
        self.assertEqual(params["certified_interval"][0]["text"], "1/1")
        self.assertEqual(params["certified_interval"][1]["text"], "11/1")

    def test_exact_derivative_rows_and_edges(self):
        payload = self.payload
        cert = payload["exact_derivative_certificate"]
        self.assertIn("N(x)=1/100-(9/5)*x^(4/5)", cert["numerator_formula"])
        self.assertEqual(cert["rational_upper_bound"]["text"], "-179/100")
        self.assertTrue(cert["upper_bound_strictly_negative"])
        self.assertIn("no floating derivative sign decision", cert["conclusion"])

        rows = payload["grid_exact_bound_rows"]
        self.assertEqual(len(rows), 11)
        self.assertTrue(all(row["strictly_negative_bound"] for row in rows))
        self.assertTrue(all(row["denominator_positive"] for row in rows))
        self.assertTrue(all(row["D_prime_certified_negative"] for row in rows))

        edges = payload["integer_edge_mean_value_rows"]
        self.assertEqual(len(edges), 10)
        self.assertTrue(all(row["therefore_D_left_greater_than_D_right"] for row in edges))
        self.assertTrue(all(not row["damping_sign_flip_possible_on_edge"] for row in edges))

    def test_summary_proof_blockers_and_hard_limits(self):
        payload = self.payload
        summary = payload["exact_rational_damping_summary"]
        self.assertTrue(summary["continuous_positive_certificate"])
        self.assertTrue(summary["continuous_strictly_decreasing_certificate"])
        self.assertEqual(summary["derivative_upper_bound_text"], "-179/100")
        self.assertTrue(summary["all_grid_bounds_negative"])
        self.assertTrue(summary["all_integer_edges_drop_by_mean_value_theorem"])
        self.assertFalse(summary["damping_can_supply_phase_sign_flips"])
        self.assertTrue(summary["phase_flip_edges_remain_phase_only"])
        self.assertTrue(summary["z2_phase_flip_edges_remain_unchanged"])
        self.assertTrue(summary["matches_previous_float_monotonicity_certificate"])

        blocker = "\n".join(payload["blocker_context"]["still_open"])
        self.assertIn("strict_damping_formula_derivation_from_nadsoliton_dynamics", blocker)
        self.assertIn("strict_phase_frequency_derivation_from_nadsoliton_dynamics", blocker)
        self.assertIn("strict_transport_derivation_from_nadsoliton_dynamics", blocker)

        proof = payload["proof_certificate"]
        self.assertIn("beta=1, eta=9/5, beta_tors=1/100", proof["parameter_step"])
        self.assertIn("-179/100<0", proof["rational_bound_step"])
        self.assertIn("mean value theorem", proof["mean_value_step"])
        self.assertIn("does not derive beta", proof["theoretical_limit"])

        hard_limits = "\n".join(payload["hard_limits"])
        self.assertIn("K_strict_gate remains the current live/full", hard_limits)
        self.assertIn("No unqualified identity K_legacy_ont == K_strict_gate", hard_limits)
        self.assertIn("No proof derives beta, eta, beta_tors, or D(x)", hard_limits)
        self.assertIn("No QW-2191 selector discharge", hard_limits)
        self.assertIn("No ToE closure", hard_limits)
        self.assertTrue(REPORT_MD.exists())


if __name__ == "__main__":
    unittest.main()
