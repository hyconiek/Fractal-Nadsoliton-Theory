from __future__ import annotations

import json
import subprocess
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "fundamental_action_reconstruction" / "scratch" / "bridge_strict_completion_phase_sign_component_quotient_exact_sequence_certificate_probe.py"
REPORT_JSON = ROOT / "fundamental_action_reconstruction" / "scratch" / "bridge_strict_completion_phase_sign_component_quotient_exact_sequence_certificate_report.json"
REPORT_MD = ROOT / "fundamental_action_reconstruction" / "scratch" / "bridge_strict_completion_phase_sign_component_quotient_exact_sequence_certificate_report.md"


class StrictCompletionPhaseSignComponentQuotientExactSequenceCertificateProbeTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        subprocess.run([sys.executable, str(SCRIPT)], check=True, cwd=ROOT)
        cls.payload = json.loads(REPORT_JSON.read_text(encoding="utf-8"))

    def test_result_kind_sources_and_exact_sequence_definition(self):
        payload = self.payload
        self.assertEqual(
            payload["result_kind"],
            "SCRATCH_STRICT_COMPLETION_PHASE_SIGN_COMPONENT_QUOTIENT_EXACT_SEQUENCE_CERTIFICATE__IMAGE_EQUALS_INTERIOR_KERNEL",
        )
        self.assertIn("imS-equals-ker", payload["status"])
        self.assertEqual(set(payload["source_reports"]), {
            "phase_sign_z2_coboundary_certificate",
            "phase_sign_component_quotient_projection_certificate",
            "phase_sign_component_quotient_lift_matrix_certificate",
        })
        self.assertIn("exact sequence", payload["grep_disambiguation"]["searched_terms"])
        definition = payload["exact_sequence_definition"]
        self.assertEqual(definition["field"], "GF(2)")
        self.assertIn("GF(2)^5 --S--> GF(2)^12", definition["sequence"])
        self.assertEqual(len(definition["F_interior_coboundary_matrix_equals_H_times_B_path"]), 7)
        self.assertEqual(len(definition["F_interior_coboundary_matrix_equals_H_times_B_path"][0]), 12)
        self.assertEqual(definition["F_times_S"], [[0, 0, 0, 0, 0]] * 7)
        self.assertEqual(len(definition["kernel_basis_rows_for_F"]), 5)
        self.assertEqual(len(definition["kernel_basis_rows_for_F"][0]), 12)

    def test_rank_nullity_reconstruction_and_exactness(self):
        payload = self.payload
        rank = payload["rank_nullity_certificate"]
        self.assertEqual(rank["rank_F_interior_coboundary"], 7)
        self.assertEqual(rank["nullity_F_interior_coboundary"], 5)
        self.assertEqual(rank["rank_S_component_expansion"], 5)
        self.assertEqual(rank["rank_kernel_basis"], 5)
        self.assertEqual(rank["rank_combined_S_columns_and_kernel_basis"], 5)
        self.assertTrue(rank["F_has_full_row_rank"])
        self.assertTrue(rank["rank_nullity_matches_component_count"])
        self.assertTrue(rank["S_rank_matches_kernel_dimension"])
        self.assertTrue(rank["combined_span_has_no_extra_dimension"])

        reconstruction = payload["reconstruction_certificate"]
        self.assertEqual(reconstruction["F_node_bits"], [0, 0, 0, 0, 0, 0, 0])
        self.assertEqual(reconstruction["Q_node_bits"], [0, 1, 0, 1, 0])
        self.assertEqual(reconstruction["S_Q_node_bits"], [0, 0, 1, 1, 1, 1, 0, 0, 1, 1, 0, 0])
        self.assertEqual(reconstruction["S_Q_node_bits"], reconstruction["projection_report_node_bits_from_S_Q_node_bits"])
        self.assertEqual(reconstruction["S_Q_node_bits"], reconstruction["lift_report_node_bits_from_S_component_bits"])

    def test_summary_proof_blockers_and_hard_limits(self):
        payload = self.payload
        summary = payload["component_quotient_exact_sequence_summary"]
        self.assertTrue(summary["F_equals_H_times_B_path_rows_by_nodes"])
        self.assertTrue(summary["F_times_S_is_zero"])
        self.assertTrue(summary["rank_F_is_interior_edge_count"])
        self.assertTrue(summary["nullity_F_is_component_count"])
        self.assertTrue(summary["S_columns_span_kernel_F"])
        self.assertTrue(summary["S_Q_projects_kernel_basis_to_itself"])
        self.assertTrue(summary["audited_node_bits_lie_in_kernel_F"])
        self.assertTrue(summary["Q_then_S_recovers_audited_node_bits"])
        self.assertTrue(summary["Q_recovers_expected_component_bits"])
        self.assertTrue(summary["matches_projection_report_round_trip"])
        self.assertTrue(summary["matches_lift_report_node_bits"])

        blocker = "\n".join(payload["blocker_context"]["still_open"])
        self.assertIn("strict_phase_frequency_derivation_from_nadsoliton_dynamics", blocker)
        self.assertIn("QW-2191_selector_obstruction", blocker)

        proof = payload["proof_certificate"]
        self.assertIn("im(S)=ker(F)", proof["exactness_step"])
        self.assertIn("S*Q fixes", proof["projector_step"])
        self.assertIn("F*b=0", proof["audited_vector_step"])
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
