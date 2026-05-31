from __future__ import annotations

import json
import subprocess
import sys
import unittest
from math import comb
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "fundamental_action_reconstruction" / "scratch" / "bridge_strict_completion_phase_sign_component_quotient_coordinate_residual_syndrome_certificate_probe.py"
REPORT_JSON = ROOT / "fundamental_action_reconstruction" / "scratch" / "bridge_strict_completion_phase_sign_component_quotient_coordinate_residual_syndrome_certificate_report.json"
REPORT_MD = ROOT / "fundamental_action_reconstruction" / "scratch" / "bridge_strict_completion_phase_sign_component_quotient_coordinate_residual_syndrome_certificate_report.md"


class StrictCompletionPhaseSignComponentQuotientCoordinateResidualSyndromeCertificateProbeTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        subprocess.run([sys.executable, str(SCRIPT)], check=True, cwd=ROOT)
        cls.payload = json.loads(REPORT_JSON.read_text(encoding="utf-8"))

    def test_result_kind_sources_and_definition(self):
        payload = self.payload
        self.assertEqual(
            payload["result_kind"],
            "SCRATCH_STRICT_COMPLETION_PHASE_SIGN_COMPONENT_QUOTIENT_COORDINATE_RESIDUAL_SYNDROME_CERTIFICATE__FULL_RESIDUAL_ENUMERATION",
        )
        self.assertIn("residual-syndrome", payload["status"])
        self.assertEqual(set(payload["source_reports"]), {
            "phase_sign_component_quotient_coordinate_isomorphism_certificate",
            "phase_sign_component_quotient_coordinate_support_minimality_certificate",
        })
        self.assertIn("coordinate residual syndrome", payload["grep_disambiguation"]["searched_terms"])
        definition = payload["coordinate_residual_definition"]
        self.assertEqual(definition["field"], "GF(2)")
        self.assertEqual(definition["coordinate_count"], 12)
        self.assertEqual(definition["node_count"], 12)
        self.assertEqual(definition["expected_coordinate_vector"], [0, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0])
        self.assertEqual(definition["zero_node_residual"], [0] * 12)

    def test_residual_syndrome_certificate(self):
        payload = self.payload
        certificate = payload["residual_syndrome_certificate"]
        self.assertEqual(certificate["coordinate_space_size"], 4096)
        self.assertEqual(certificate["unique_residual_syndrome_count"], 4096)
        self.assertEqual(certificate["minimum_nonzero_residual_weight"], 1)
        self.assertEqual(len(certificate["zero_syndrome_rows"]), 1)
        zero_row = certificate["zero_syndrome_rows"][0]
        self.assertEqual(zero_row["coordinate_vector"], [0, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0])
        self.assertEqual(zero_row["residual_syndrome"], [0] * 12)
        self.assertEqual(zero_row["residual_weight"], 0)
        self.assertTrue(certificate["minimum_nonzero_residual_examples"])
        self.assertTrue(certificate["nonzero_failure_examples"])
        expected_histogram = {str(weight): comb(12, weight) for weight in range(13)}
        self.assertEqual(certificate["residual_weight_histogram"], expected_histogram)
        self.assertEqual(certificate["expected_binomial_residual_weight_histogram"], expected_histogram)

    def test_summary_proof_blockers_and_hard_limits(self):
        payload = self.payload
        summary = payload["coordinate_residual_syndrome_summary"]
        self.assertTrue(summary["enumerated_all_2^12_coordinate_vectors"])
        self.assertTrue(summary["all_residual_syndromes_unique"])
        self.assertTrue(summary["zero_syndrome_unique"])
        self.assertTrue(summary["zero_syndrome_coordinate_is_audited"])
        self.assertTrue(summary["every_nonmatching_coordinate_has_nonzero_residual"])
        self.assertTrue(summary["minimum_nonzero_residual_weight_is_1"])
        self.assertTrue(summary["residual_weight_histogram_is_binomial"])
        self.assertTrue(summary["coordinate_delta_histogram_sums_to_coordinate_space"])
        self.assertTrue(summary["matches_coordinate_support_minimality_unique_vector"])

        blocker = "\n".join(payload["blocker_context"]["still_open"])
        self.assertIn("strict_phase_frequency_derivation_from_nadsoliton_dynamics", blocker)
        self.assertIn("QW-2191_selector_obstruction", blocker)

        proof = payload["proof_certificate"]
        self.assertIn("residual_syndrome(c)=U*c+target_node_bits", proof["residual_step"])
        self.assertIn("Exactly one row", proof["zero_step"])
        self.assertIn("4095 coordinate vectors", proof["failure_step"])
        self.assertIn("C(12,k)", proof["histogram_step"])
        self.assertIn("does not derive omega/phi", proof["theoretical_limit"])

        hard_limits = "\n".join(payload["hard_limits"])
        self.assertIn("K_strict_gate remains the current live/full", hard_limits)
        self.assertIn("No unqualified identity K_legacy_ont == K_strict_gate", hard_limits)
        self.assertIn("No beta_tors -> chi_11 theorem", hard_limits)
        self.assertIn("No QW-2191 selector discharge", hard_limits)
        self.assertIn("No ToE closure", hard_limits)
        self.assertTrue(REPORT_MD.exists())


if __name__ == "__main__":
    unittest.main()
