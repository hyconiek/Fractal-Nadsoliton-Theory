from __future__ import annotations

import json
import subprocess
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "fundamental_action_reconstruction" / "scratch" / "bridge_strict_completion_phase_sign_component_quotient_complement_inverse_certificate_probe.py"
REPORT_JSON = ROOT / "fundamental_action_reconstruction" / "scratch" / "bridge_strict_completion_phase_sign_component_quotient_complement_inverse_certificate_report.json"
REPORT_MD = ROOT / "fundamental_action_reconstruction" / "scratch" / "bridge_strict_completion_phase_sign_component_quotient_complement_inverse_certificate_report.md"


class StrictCompletionPhaseSignComponentQuotientComplementInverseCertificateProbeTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        subprocess.run([sys.executable, str(SCRIPT)], check=True, cwd=ROOT)
        cls.payload = json.loads(REPORT_JSON.read_text(encoding="utf-8"))

    def test_result_kind_sources_and_complement_definition(self):
        payload = self.payload
        self.assertEqual(
            payload["result_kind"],
            "SCRATCH_STRICT_COMPLETION_PHASE_SIGN_COMPONENT_QUOTIENT_COMPLEMENT_INVERSE_CERTIFICATE__DIRECT_SUM_AND_INTERIOR_INVERSE",
        )
        self.assertIn("complement-direct-sum", payload["status"])
        self.assertEqual(set(payload["source_reports"]), {
            "phase_sign_component_quotient_projection_certificate",
            "phase_sign_component_quotient_exact_sequence_certificate",
            "phase_sign_z2_coboundary_certificate",
        })
        self.assertIn("component quotient complement", payload["grep_disambiguation"]["searched_terms"])
        definition = payload["complement_definition"]
        self.assertEqual(definition["field"], "GF(2)")
        self.assertEqual(len(definition["N_residual_basis_matrix_nodes_by_residuals"]), 12)
        self.assertEqual(len(definition["N_residual_basis_matrix_nodes_by_residuals"][0]), 7)
        self.assertEqual(len(definition["residual_basis_rows"]), 7)
        self.assertEqual(definition["Q_times_N"], [[0] * 7] * 5)
        self.assertEqual(definition["F_times_S"], [[0] * 5] * 7)
        self.assertEqual(definition["FN_inverse_times_FN"], [[1 if row == col else 0 for col in range(7)] for row in range(7)])
        self.assertEqual(definition["FN_times_FN_inverse"], [[1 if row == col else 0 for col in range(7)] for row in range(7)])

    def test_rank_reconstruction_and_direct_sum(self):
        payload = self.payload
        rank = payload["rank_certificate"]
        self.assertEqual(rank["rank_S"], 5)
        self.assertEqual(rank["rank_N"], 7)
        self.assertEqual(rank["rank_block_basis_S_then_N"], 12)
        self.assertEqual(rank["rank_F_times_N"], 7)
        self.assertTrue(rank["S_and_N_span_all_nodes"])
        self.assertTrue(rank["N_has_residual_rank"])
        self.assertTrue(rank["F_times_N_invertible"])

        reconstruction = payload["reconstruction_certificate"]
        self.assertEqual(reconstruction["component_bits_Q_node_bits"], [0, 1, 0, 1, 0])
        self.assertEqual(reconstruction["quotient_part_S_Q_node_bits"], [0, 0, 1, 1, 1, 1, 0, 0, 1, 1, 0, 0])
        self.assertEqual(reconstruction["residual_part_node_bits_plus_SQnode_bits"], [0] * 12)
        self.assertEqual(reconstruction["F_node_bits"], [0] * 7)
        self.assertEqual(reconstruction["residual_coordinates_from_FN_inverse_F_node_bits"], [0] * 7)
        self.assertEqual(reconstruction["residual_from_coordinates"], [0] * 12)

    def test_summary_proof_blockers_and_hard_limits(self):
        payload = self.payload
        summary = payload["component_quotient_complement_inverse_summary"]
        self.assertTrue(summary["Q_times_S_is_identity"])
        self.assertTrue(summary["Q_times_N_is_zero"])
        self.assertTrue(summary["F_times_S_is_zero"])
        self.assertTrue(summary["F_times_N_is_invertible"])
        self.assertTrue(summary["block_basis_is_full_rank_12"])
        self.assertTrue(summary["direct_sum_dimensions_add_to_12"])
        self.assertTrue(summary["audited_node_bits_have_zero_residual_part"])
        self.assertTrue(summary["audited_residual_coordinates_zero"])
        self.assertTrue(summary["audited_quotient_part_matches_node_bits"])
        self.assertTrue(summary["component_bits_match_expected"])
        self.assertTrue(summary["matches_exact_sequence_F_node_bits"])

        blocker = "\n".join(payload["blocker_context"]["still_open"])
        self.assertIn("strict_phase_frequency_derivation_from_nadsoliton_dynamics", blocker)
        self.assertIn("QW-2191_selector_obstruction", blocker)

        proof = payload["proof_certificate"]
        self.assertIn("[S N] has rank 12", proof["direct_sum_step"])
        self.assertIn("Q*N=0", proof["annihilation_step"])
        self.assertIn("F*N has rank 7", proof["inverse_step"])
        self.assertIn("residual coordinates are all zero", proof["audited_vector_step"])
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
