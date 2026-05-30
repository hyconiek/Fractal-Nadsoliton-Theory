from __future__ import annotations

import json
import subprocess
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "fundamental_action_reconstruction" / "scratch" / "bridge_phase_normalized_z12_inverse_branch_interval_certificate_probe.py"
REPORT_JSON = ROOT / "fundamental_action_reconstruction" / "scratch" / "bridge_phase_normalized_z12_inverse_branch_interval_certificate_report.json"
REPORT_MD = ROOT / "fundamental_action_reconstruction" / "scratch" / "bridge_phase_normalized_z12_inverse_branch_interval_certificate_report.md"


class PhaseNormalizedZ12InverseBranchIntervalCertificateProbeTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        subprocess.run([sys.executable, str(SCRIPT)], check=True, cwd=ROOT)
        cls.payload = json.loads(REPORT_JSON.read_text(encoding="utf-8"))

    def test_result_kind_constants_and_interval_certificate(self):
        payload = self.payload
        self.assertEqual(
            payload["result_kind"],
            "SCRATCH_PHASE_NORMALIZED_Z12_INVERSE_BRANCH_INTERVAL_CERTIFICATE__OUTPUT_MATCHING_NOT_STRICT_DERIVATION",
        )
        self.assertEqual(payload["status"], "finite-z12-monotone-inverse-branch-certified-by-interval-brackets")
        constants = payload["constants"]
        self.assertEqual(constants["branch_interval"], [0.0, 2.0])
        self.assertEqual(constants["domain_d_values"], list(range(12)))
        self.assertEqual(constants["bisection_steps"], 90)
        self.assertEqual(constants["target_q_power"], "256/243")
        self.assertEqual(constants["target_eta"], "9/5")

        cert = payload["interval_certificate"]
        self.assertTrue(cert["all_targets_inside_legacy_interval_range"])
        self.assertLess(cert["legacy_derivative_upper_bound_on_interval"], 0.0)
        self.assertTrue(cert["legacy_derivative_bound_strictly_negative"])
        self.assertTrue(cert["all_root_intervals_sign_bracket_or_endpoint"])
        self.assertLess(cert["max_root_interval_width"], 1e-12)
        self.assertLess(cert["max_mid_residual_abs"], 1e-12)
        self.assertTrue(cert["strict_targets_strictly_decreasing_on_z12"])
        self.assertGreater(cert["min_strict_target_drop_margin"], 0.0)
        self.assertTrue(cert["roots_monotone_non_decreasing_on_z12"])
        self.assertGreater(cert["min_root_increment"], 0.0)
        self.assertTrue(cert["all_midpoints_in_first_legacy_interval"])
        self.assertAlmostEqual(cert["z12_to_legacy_compression_x11_over_11"], 0.12201655847945313)

    def test_root_rows_are_bracketed_and_monotone(self):
        rows = self.payload["root_rows"]
        self.assertEqual(len(rows), 12)
        roots = [row["root_interval"]["mid"] for row in rows]
        targets = [row["strict_norm_target"] for row in rows]
        self.assertEqual(rows[0]["d"], 0)
        self.assertEqual(rows[-1]["d"], 11)
        self.assertAlmostEqual(roots[0], 0.0)
        self.assertTrue(all(left < right for left, right in zip(roots, roots[1:])))
        self.assertTrue(all(left > right for left, right in zip(targets, targets[1:])))
        for row in rows:
            self.assertTrue(row["mid_in_first_legacy_interval"])
            self.assertTrue(row["root_interval"]["sign_bracket_or_endpoint"])
            self.assertLessEqual(row["root_interval"]["width"], 1e-12)
            self.assertAlmostEqual(row["legacy_norm_at_mid"], row["strict_norm_target"], delta=1e-12)

    def test_upstream_update_proof_guardrails_and_markdown(self):
        update = self.payload["upstream_blocker_update"]
        self.assertIn("OPEN_MONOTONE_INVERSE_BRANCH", update["previous_monotone_status"])
        self.assertFalse(update["previous_global_admissibility_note"])
        self.assertEqual(update["blocker_lattice_global_z12_status_before_this_probe"], "OPEN_BLOCKER")
        self.assertIn("global_z12_output_matching_map_certified", update["refined_status_after_this_probe"])
        self.assertIn("strict_transport_derivation", update["refined_status_after_this_probe"])
        self.assertIn("orientation_chi11_source", update["remaining_bridge_blockers"])

        proof = self.payload["exact_proof_certificate"]
        self.assertIn("IVT", proof["existence"])
        self.assertIn("derivative upper bound", proof["uniqueness"])
        self.assertIn("finite strict targets", proof["monotonicity"])
        self.assertIn("does not derive", proof["blocker_refinement"])

        ontology = self.payload["ontology_guardrail"]
        self.assertIn("nadsoliton itself", ontology["allowed_reading"])
        self.assertIn("No separate informational layer", ontology["forbidden_reading"])

        hard_limits = "\n".join(self.payload["hard_limits"])
        self.assertIn("No unqualified identity K_legacy_ont == K_strict_gate", hard_limits)
        self.assertIn("No strict-side derivation", hard_limits)
        self.assertIn("No beta_tors -> chi_11 theorem", hard_limits)
        self.assertIn("No legacy physical-role transfer", hard_limits)
        self.assertIn("No QW-2191 discharge", hard_limits)
        self.assertIn("No ToE closure", hard_limits)
        self.assertTrue(REPORT_MD.exists())


if __name__ == "__main__":
    unittest.main()
