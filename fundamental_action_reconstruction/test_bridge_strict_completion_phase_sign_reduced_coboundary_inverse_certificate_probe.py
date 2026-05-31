from __future__ import annotations

import json
import subprocess
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "fundamental_action_reconstruction" / "scratch" / "bridge_strict_completion_phase_sign_reduced_coboundary_inverse_certificate_probe.py"
REPORT_JSON = ROOT / "fundamental_action_reconstruction" / "scratch" / "bridge_strict_completion_phase_sign_reduced_coboundary_inverse_certificate_report.json"
REPORT_MD = ROOT / "fundamental_action_reconstruction" / "scratch" / "bridge_strict_completion_phase_sign_reduced_coboundary_inverse_certificate_report.md"


class StrictCompletionPhaseSignReducedCoboundaryInverseCertificateProbeTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        subprocess.run([sys.executable, str(SCRIPT)], check=True, cwd=ROOT)
        cls.payload = json.loads(REPORT_JSON.read_text(encoding="utf-8"))

    def test_result_kind_sources_and_reduced_complex_definition(self):
        payload = self.payload
        self.assertEqual(
            payload["result_kind"],
            "SCRATCH_STRICT_COMPLETION_PHASE_SIGN_REDUCED_COBOUNDARY_INVERSE_CERTIFICATE__ANCHORED_GF2_ISOMORPHISM",
        )
        self.assertIn("anchored-reduced-coboundary", payload["status"])
        self.assertEqual(set(payload["source_reports"]), {
            "phase_sign_z2_coboundary_certificate",
            "phase_sign_gf2_linear_system_certificate",
            "phase_sign_path_cohomology_triviality_certificate",
            "phase_sign_cycle_closure_obstruction_certificate",
        })
        self.assertIn("reduced coboundary", payload["grep_disambiguation"]["searched_terms"])
        definition = payload["anchored_reduced_complex_definition"]
        self.assertEqual(definition["field"], "GF(2)")
        self.assertEqual(definition["anchor_node"], 0)
        self.assertEqual(len(definition["reduced_coboundary_R_matrix_anchor_column_removed"]), 11)
        self.assertEqual(len(definition["prefix_inverse_P_matrix"]), 11)
        self.assertEqual(definition["reduced_coboundary_R_matrix_anchor_column_removed"][0], [1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0])
        self.assertEqual(definition["prefix_inverse_P_matrix"][10], [1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1])

    def test_rank_inverse_reconstruction_and_rows(self):
        payload = self.payload
        rank = payload["rank_certificate"]
        self.assertEqual(rank["rank"], 11)
        self.assertEqual(rank["nullity"], 0)
        self.assertEqual(rank["determinant_mod2"], 1)
        self.assertTrue(rank["left_block_is_identity"])
        self.assertEqual([row["pivot_col"] for row in rank["pivot_rows"]], list(range(11)))

        inverse = payload["two_sided_inverse_certificate"]
        self.assertTrue(inverse["R_times_P_is_identity"])
        self.assertTrue(inverse["P_times_R_is_identity"])
        self.assertEqual(inverse["R_times_P"], inverse["identity_11"])
        self.assertEqual(inverse["P_times_R"], inverse["identity_11"])

        reconstruction = payload["reconstruction_certificate"]
        self.assertEqual(reconstruction["full_node_bits_from_anchor_and_reduced_inverse"], [0, 0, 1, 1, 1, 1, 0, 0, 1, 1, 0, 0])
        self.assertEqual(reconstruction["edge_bits_from_R_tail_node_bits"], [0, 1, 0, 0, 0, 1, 0, 1, 0, 1, 0])
        self.assertEqual(reconstruction["flip_edges_from_edge_bits"], ["1->2", "5->6", "7->8", "9->10"])

        self.assertEqual(len(payload["reduced_equation_rows"]), 11)
        self.assertTrue(all(row["matches_audited_edge_bit"] for row in payload["reduced_equation_rows"]))
        self.assertTrue(payload["reduced_equation_rows"][0]["anchor_column_removed"])
        self.assertEqual(len(payload["prefix_inverse_rows"]), 11)
        self.assertTrue(all(row["matches_audited_tail_node_bit"] for row in payload["prefix_inverse_rows"]))

    def test_summary_proof_blockers_and_hard_limits(self):
        payload = self.payload
        summary = payload["reduced_coboundary_inverse_summary"]
        self.assertEqual(summary["rank_R"], 11)
        self.assertEqual(summary["nullity_R"], 0)
        self.assertEqual(summary["determinant_mod2_R"], 1)
        self.assertTrue(summary["anchored_map_is_isomorphism"])
        self.assertTrue(summary["two_sided_prefix_inverse_verified"])
        self.assertTrue(summary["full_node_bits_match_z2"])
        self.assertTrue(summary["edge_bits_match_z2"])
        self.assertTrue(summary["matches_gf2_linear_system_prefix_inverse"])
        self.assertTrue(summary["matches_gf2_linear_system_solution"])
        self.assertTrue(summary["matches_path_cohomology_anchor_reconstruction"])
        self.assertTrue(summary["inherits_path_h1_zero"])
        self.assertTrue(summary["contrasts_cycle_h1_one_boundary_check"])
        self.assertTrue(summary["all_reduced_equations_pass"])
        self.assertTrue(summary["all_prefix_inverse_rows_pass"])

        blocker = "\n".join(payload["blocker_context"]["still_open"])
        self.assertIn("strict_phase_frequency_derivation_from_nadsoliton_dynamics", blocker)
        self.assertIn("QW-2191_selector_obstruction", blocker)

        proof = payload["proof_certificate"]
        self.assertIn("anchor column b(0)", proof["anchor_step"])
        self.assertIn("rank(R)=11", proof["rank_step"])
        self.assertIn("R*P=I_11", proof["inverse_step"])
        self.assertIn("1->2, 5->6, 7->8, and 9->10", proof["reconstruction_step"])
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
