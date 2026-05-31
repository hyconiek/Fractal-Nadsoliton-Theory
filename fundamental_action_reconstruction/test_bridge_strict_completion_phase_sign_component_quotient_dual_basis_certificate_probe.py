from __future__ import annotations

import json
import subprocess
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "fundamental_action_reconstruction" / "scratch" / "bridge_strict_completion_phase_sign_component_quotient_dual_basis_certificate_probe.py"
REPORT_JSON = ROOT / "fundamental_action_reconstruction" / "scratch" / "bridge_strict_completion_phase_sign_component_quotient_dual_basis_certificate_report.json"
REPORT_MD = ROOT / "fundamental_action_reconstruction" / "scratch" / "bridge_strict_completion_phase_sign_component_quotient_dual_basis_certificate_report.md"


class StrictCompletionPhaseSignComponentQuotientDualBasisCertificateProbeTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        subprocess.run([sys.executable, str(SCRIPT)], check=True, cwd=ROOT)
        cls.payload = json.loads(REPORT_JSON.read_text(encoding="utf-8"))

    def test_result_kind_sources_and_dual_definition(self):
        payload = self.payload
        self.assertEqual(
            payload["result_kind"],
            "SCRATCH_STRICT_COMPLETION_PHASE_SIGN_COMPONENT_QUOTIENT_DUAL_BASIS_CERTIFICATE__PAIRING_AND_DECOMPOSITION",
        )
        self.assertIn("dual-basis", payload["status"])
        self.assertEqual(set(payload["source_reports"]), {
            "phase_sign_component_quotient_coordinate_isomorphism_certificate",
            "phase_sign_z2_coboundary_certificate",
        })
        self.assertIn("dual basis", payload["grep_disambiguation"]["searched_terms"])
        definition = payload["dual_basis_definition"]
        self.assertEqual(definition["field"], "GF(2)")
        self.assertEqual(len(definition["coordinate_functional_rows_T_i"]), 12)
        self.assertEqual(len(definition["basis_columns_U_i"]), 12)
        self.assertEqual(definition["pairing_matrix_T_i_on_U_j"], [[1 if row == col else 0 for col in range(12)] for row in range(12)])
        self.assertEqual(definition["active_coordinate_labels"], ["quotient_component_1", "quotient_component_3"])
        self.assertEqual(len(definition["decomposition_rows"]), 12)

    def test_rank_pairing_and_reconstruction(self):
        payload = self.payload
        rank = payload["rank_pairing_certificate"]
        self.assertEqual(rank["rank_T_functionals"], 12)
        self.assertEqual(rank["rank_U_basis"], 12)
        self.assertEqual(rank["rank_pairing_matrix"], 12)
        self.assertTrue(rank["pairing_matrix_is_identity"])
        self.assertTrue(rank["T_functionals_full_rank"])
        self.assertTrue(rank["U_basis_full_rank"])

        reconstruction = payload["reconstruction_certificate"]
        self.assertEqual(reconstruction["coordinate_vector_T_b"], [0, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0])
        self.assertEqual(reconstruction["reconstructed_from_dual_basis"], [0, 0, 1, 1, 1, 1, 0, 0, 1, 1, 0, 0])
        self.assertEqual(reconstruction["coordinate_vector_T_b"], reconstruction["coordinate_report_coordinate_vector"])
        self.assertEqual(reconstruction["reconstructed_from_dual_basis"], reconstruction["coordinate_report_reconstructed_node_bits"])

    def test_summary_proof_blockers_and_hard_limits(self):
        payload = self.payload
        summary = payload["component_quotient_dual_basis_summary"]
        self.assertTrue(summary["pairing_matrix_is_identity"])
        self.assertTrue(summary["T_and_U_full_rank"])
        self.assertTrue(summary["coordinate_vector_matches_expected"])
        self.assertTrue(summary["active_coordinates_are_quotient_components_1_and_3"])
        self.assertTrue(summary["all_interior_residual_coordinates_zero"])
        self.assertTrue(summary["dual_basis_reconstructs_node_bits"])
        self.assertTrue(summary["matches_coordinate_isomorphism_report"])

        blocker = "\n".join(payload["blocker_context"]["still_open"])
        self.assertIn("strict_phase_frequency_derivation_from_nadsoliton_dynamics", blocker)
        self.assertIn("QW-2191_selector_obstruction", blocker)

        proof = payload["proof_certificate"]
        self.assertIn("<T_i,U_j>=delta_ij", proof["pairing_step"])
        self.assertIn("quotient_component_1 and quotient_component_3", proof["decomposition_step"])
        self.assertIn("seven interior residual coordinates vanish", proof["residual_step"])
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
