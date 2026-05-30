from __future__ import annotations

import json
import subprocess
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "fundamental_action_reconstruction" / "scratch" / "bridge_strict_completion_phase_sign_gf2_linear_system_certificate_probe.py"
REPORT_JSON = ROOT / "fundamental_action_reconstruction" / "scratch" / "bridge_strict_completion_phase_sign_gf2_linear_system_certificate_report.json"
REPORT_MD = ROOT / "fundamental_action_reconstruction" / "scratch" / "bridge_strict_completion_phase_sign_gf2_linear_system_certificate_report.md"


class StrictCompletionPhaseSignGF2LinearSystemCertificateProbeTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        subprocess.run([sys.executable, str(SCRIPT)], check=True, cwd=ROOT)
        cls.payload = json.loads(REPORT_JSON.read_text(encoding="utf-8"))

    def test_result_kind_sources_and_system_definition(self):
        payload = self.payload
        self.assertEqual(
            payload["result_kind"],
            "SCRATCH_STRICT_COMPLETION_PHASE_SIGN_GF2_LINEAR_SYSTEM_CERTIFICATE__FULL_RANK_PREFIX_RECONSTRUCTION",
        )
        self.assertIn("full-rank", payload["status"])
        self.assertEqual(set(payload["source_reports"]), {
            "phase_sign_z2_coboundary_certificate",
            "phase_sign_edge_support_minimality_certificate",
        })
        self.assertIn("phase_sign_gf2", payload["grep_disambiguation"]["searched_terms"])
        system = payload["linear_system_definition"]
        self.assertEqual(system["field"], "GF(2)")
        self.assertEqual(system["anchor_bit_b0"], 0)
        self.assertEqual(len(system["prefix_matrix_L"]), 11)
        self.assertEqual(system["prefix_matrix_L"][0], [1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0])
        self.assertEqual(system["prefix_matrix_L"][10], [1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1])
        self.assertIn("first-difference", system["explicit_inverse_first_difference_matrix"] and payload["proof_certificate"]["inverse_step"])

    def test_row_reduction_inverse_equations_and_solution(self):
        payload = self.payload
        reduction = payload["gf2_row_reduction_certificate"]
        self.assertEqual(reduction["rank"], 11)
        self.assertEqual(reduction["nullity"], 0)
        self.assertTrue(reduction["left_block_is_identity"])
        self.assertEqual([row["pivot_col"] for row in reduction["pivot_rows"]], list(range(11)))
        self.assertEqual(reduction["solution"], [0, 1, 0, 0, 0, 1, 0, 1, 0, 1, 0])

        inverse = payload["inverse_checks"]
        self.assertTrue(inverse["L_times_inverse_is_identity"])
        self.assertTrue(inverse["inverse_times_L_is_identity"])
        self.assertEqual(inverse["inverse_solution_edge_bits"], reduction["solution"])
        self.assertTrue(inverse["inverse_solution_matches_row_reduction_solution"])

        equation_rows = payload["equation_check_rows"]
        self.assertEqual(len(equation_rows), 11)
        self.assertTrue(all(row["equation_passes"] for row in equation_rows))
        self.assertEqual(equation_rows[1]["prefix_edges"], ["0->1", "1->2"])

        solution_rows = payload["solution_edge_rows"]
        self.assertEqual(len(solution_rows), 11)
        self.assertTrue(all(row["matches_z2_edge_bit"] for row in solution_rows))
        self.assertEqual([row["edge"] for row in solution_rows if row["is_flip_edge"]], ["1->2", "5->6", "7->8", "9->10"])

    def test_summary_proof_blockers_and_hard_limits(self):
        payload = self.payload
        summary = payload["gf2_linear_system_summary"]
        self.assertEqual(summary["rank"], 11)
        self.assertEqual(summary["nullity"], 0)
        self.assertEqual(summary["determinant_mod2"], 1)
        self.assertTrue(summary["unique_solution"])
        self.assertEqual(summary["solution_hamming_weight"], 4)
        self.assertEqual(summary["solution_flip_edges"], ["1->2", "5->6", "7->8", "9->10"])
        self.assertTrue(summary["all_equations_pass"])
        self.assertTrue(summary["all_solution_edges_match_z2"])
        self.assertTrue(summary["matches_expected_node_bits"])
        self.assertTrue(summary["matches_expected_edge_bits"])
        self.assertTrue(summary["matches_expected_flip_edges"])
        self.assertTrue(summary["matches_edge_support_minimality_solution"])
        self.assertTrue(summary["inherits_edge_support_minimality"])

        blocker = "\n".join(payload["blocker_context"]["still_open"])
        self.assertIn("strict_phase_frequency_derivation_from_nadsoliton_dynamics", blocker)
        self.assertIn("QW-2191_selector_obstruction", blocker)

        proof = payload["proof_certificate"]
        self.assertIn("L e = b_tail xor b0", proof["system_step"])
        self.assertIn("rank(L)=11", proof["rank_step"])
        self.assertIn("two-sided inverse", proof["inverse_step"])
        self.assertIn("1->2, 5->6, 7->8, and 9->10", proof["solution_step"])
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
