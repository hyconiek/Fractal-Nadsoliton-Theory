from __future__ import annotations

import json
import subprocess
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "fundamental_action_reconstruction" / "scratch" / "bridge_strict_completion_phase_sign_component_quotient_coordinate_support_minimality_certificate_probe.py"
REPORT_JSON = ROOT / "fundamental_action_reconstruction" / "scratch" / "bridge_strict_completion_phase_sign_component_quotient_coordinate_support_minimality_certificate_report.json"
REPORT_MD = ROOT / "fundamental_action_reconstruction" / "scratch" / "bridge_strict_completion_phase_sign_component_quotient_coordinate_support_minimality_certificate_report.md"


class StrictCompletionPhaseSignComponentQuotientCoordinateSupportMinimalityCertificateProbeTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        subprocess.run([sys.executable, str(SCRIPT)], check=True, cwd=ROOT)
        cls.payload = json.loads(REPORT_JSON.read_text(encoding="utf-8"))

    def test_result_kind_sources_and_definition(self):
        payload = self.payload
        self.assertEqual(
            payload["result_kind"],
            "SCRATCH_STRICT_COMPLETION_PHASE_SIGN_COMPONENT_QUOTIENT_COORDINATE_SUPPORT_MINIMALITY_CERTIFICATE__FULL_COORDINATE_ENUMERATION",
        )
        self.assertIn("support-minimality", payload["status"])
        self.assertEqual(set(payload["source_reports"]), {
            "phase_sign_component_quotient_dual_basis_certificate",
            "phase_sign_component_quotient_coordinate_isomorphism_certificate",
        })
        self.assertIn("coordinate support minimality", payload["grep_disambiguation"]["searched_terms"])
        definition = payload["coordinate_support_definition"]
        self.assertEqual(definition["field"], "GF(2)")
        self.assertEqual(definition["coordinate_count"], 12)
        self.assertEqual(definition["expected_coordinate_vector"], [0, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0])
        self.assertEqual(definition["expected_active_coordinate_labels"], ["quotient_component_1", "quotient_component_3"])

    def test_enumeration_certificate(self):
        payload = self.payload
        enumeration = payload["enumeration_certificate"]
        self.assertEqual(enumeration["coordinate_space_size"], 4096)
        self.assertEqual(sum(enumeration["weight_histogram"].values()), 4096)
        self.assertEqual(enumeration["min_matching_weight"], 2)
        self.assertEqual(len(enumeration["matching_rows"]), 1)
        row = enumeration["matching_rows"][0]
        self.assertEqual(row["coordinate_vector"], [0, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0])
        self.assertEqual(row["support_weight"], 2)
        self.assertEqual(row["support_labels"], ["quotient_component_1", "quotient_component_3"])
        self.assertTrue(enumeration["lower_weight_failure_examples"])
        self.assertTrue(enumeration["same_weight_failure_examples"])

    def test_summary_proof_blockers_and_hard_limits(self):
        payload = self.payload
        summary = payload["coordinate_support_minimality_summary"]
        self.assertTrue(summary["enumerated_all_2^12_coordinate_vectors"])
        self.assertTrue(summary["unique_matching_coordinate_vector"])
        self.assertTrue(summary["matching_coordinate_vector_equals_dual_basis"])
        self.assertTrue(summary["minimum_matching_weight_is_2"])
        self.assertTrue(summary["all_lower_weight_vectors_fail"])
        self.assertTrue(summary["all_other_weight_2_vectors_fail"])
        self.assertTrue(summary["active_coordinates_are_expected_quotient_components"])
        self.assertTrue(summary["all_interior_residual_coordinates_zero"])
        self.assertTrue(summary["weight_histogram_sums_to_coordinate_space"])

        blocker = "\n".join(payload["blocker_context"]["still_open"])
        self.assertIn("strict_phase_frequency_derivation_from_nadsoliton_dynamics", blocker)
        self.assertIn("QW-2191_selector_obstruction", blocker)

        proof = payload["proof_certificate"]
        self.assertIn("All 2^12 coordinate vectors", proof["enumeration_step"])
        self.assertIn("Exactly one coordinate vector", proof["uniqueness_step"])
        self.assertIn("weight-0 or weight-1", proof["minimality_step"])
        self.assertIn("quotient_component_1 plus quotient_component_3", proof["support_step"])
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
