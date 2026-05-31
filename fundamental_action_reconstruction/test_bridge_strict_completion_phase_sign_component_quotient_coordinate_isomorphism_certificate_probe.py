from __future__ import annotations

import json
import subprocess
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "fundamental_action_reconstruction" / "scratch" / "bridge_strict_completion_phase_sign_component_quotient_coordinate_isomorphism_certificate_probe.py"
REPORT_JSON = ROOT / "fundamental_action_reconstruction" / "scratch" / "bridge_strict_completion_phase_sign_component_quotient_coordinate_isomorphism_certificate_report.json"
REPORT_MD = ROOT / "fundamental_action_reconstruction" / "scratch" / "bridge_strict_completion_phase_sign_component_quotient_coordinate_isomorphism_certificate_report.md"


class StrictCompletionPhaseSignComponentQuotientCoordinateIsomorphismCertificateProbeTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        subprocess.run([sys.executable, str(SCRIPT)], check=True, cwd=ROOT)
        cls.payload = json.loads(REPORT_JSON.read_text(encoding="utf-8"))

    def test_result_kind_sources_and_coordinate_definition(self):
        payload = self.payload
        self.assertEqual(
            payload["result_kind"],
            "SCRATCH_STRICT_COMPLETION_PHASE_SIGN_COMPONENT_QUOTIENT_COORDINATE_ISOMORPHISM_CERTIFICATE__TWO_SIDED_BLOCK_INVERSE",
        )
        self.assertIn("coordinate-isomorphism", payload["status"])
        self.assertEqual(set(payload["source_reports"]), {
            "phase_sign_component_quotient_complement_inverse_certificate",
            "phase_sign_z2_coboundary_certificate",
        })
        self.assertIn("coordinate isomorphism", payload["grep_disambiguation"]["searched_terms"])
        definition = payload["coordinate_isomorphism_definition"]
        self.assertEqual(definition["field"], "GF(2)")
        self.assertEqual(len(definition["coordinate_map_T_rows_quotient_then_interior"]), 12)
        self.assertEqual(len(definition["coordinate_map_T_rows_quotient_then_interior"][0]), 12)
        self.assertEqual(len(definition["inverse_map_U_columns_quotient_then_interior"]), 12)
        self.assertEqual(len(definition["inverse_map_U_columns_quotient_then_interior"][0]), 12)
        self.assertEqual(definition["T_times_U"], [[1 if row == col else 0 for col in range(12)] for row in range(12)])
        self.assertEqual(definition["U_times_T"], [[1 if row == col else 0 for col in range(12)] for row in range(12)])

    def test_rank_inverse_and_reconstruction(self):
        payload = self.payload
        rank = payload["rank_inverse_certificate"]
        self.assertEqual(rank["rank_T_coordinate_map"], 12)
        self.assertEqual(rank["rank_U_inverse_map"], 12)
        self.assertTrue(rank["T_has_full_rank"])
        self.assertTrue(rank["U_has_full_rank"])
        self.assertTrue(rank["T_times_U_is_identity"])
        self.assertTrue(rank["U_times_T_is_identity"])

        reconstruction = payload["reconstruction_certificate"]
        self.assertEqual(reconstruction["quotient_component_coordinates"], [0, 1, 0, 1, 0])
        self.assertEqual(reconstruction["interior_residual_coordinates"], [0] * 7)
        self.assertEqual(reconstruction["node_bits_from_U_T_node_bits"], [0, 0, 1, 1, 1, 1, 0, 0, 1, 1, 0, 0])
        self.assertEqual(reconstruction["quotient_component_coordinates"], reconstruction["complement_report_component_bits"])
        self.assertEqual(reconstruction["interior_residual_coordinates"], reconstruction["complement_report_residual_coordinates"])

    def test_summary_proof_blockers_and_hard_limits(self):
        payload = self.payload
        summary = payload["component_quotient_coordinate_isomorphism_summary"]
        self.assertTrue(summary["T_rank_full_12"])
        self.assertTrue(summary["U_rank_full_12"])
        self.assertTrue(summary["T_U_identity"])
        self.assertTrue(summary["U_T_identity"])
        self.assertTrue(summary["coordinate_vector_splits_expected"])
        self.assertTrue(summary["U_reconstructs_audited_node_bits"])
        self.assertTrue(summary["matches_complement_report_coordinates"])

        blocker = "\n".join(payload["blocker_context"]["still_open"])
        self.assertIn("strict_phase_frequency_derivation_from_nadsoliton_dynamics", blocker)
        self.assertIn("QW-2191_selector_obstruction", blocker)

        proof = payload["proof_certificate"]
        self.assertIn("T=[Q;F]", proof["coordinate_step"])
        self.assertIn("T*U=I_12", proof["inverse_step"])
        self.assertIn("0,0,0,0,0,0,0", proof["audited_vector_step"])
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
