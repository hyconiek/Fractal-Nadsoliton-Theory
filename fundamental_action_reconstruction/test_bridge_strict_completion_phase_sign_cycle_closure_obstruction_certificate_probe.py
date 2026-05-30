from __future__ import annotations

import json
import subprocess
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "fundamental_action_reconstruction" / "scratch" / "bridge_strict_completion_phase_sign_cycle_closure_obstruction_certificate_probe.py"
REPORT_JSON = ROOT / "fundamental_action_reconstruction" / "scratch" / "bridge_strict_completion_phase_sign_cycle_closure_obstruction_certificate_report.json"
REPORT_MD = ROOT / "fundamental_action_reconstruction" / "scratch" / "bridge_strict_completion_phase_sign_cycle_closure_obstruction_certificate_report.md"


class StrictCompletionPhaseSignCycleClosureObstructionCertificateProbeTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        subprocess.run([sys.executable, str(SCRIPT)], check=True, cwd=ROOT)
        cls.payload = json.loads(REPORT_JSON.read_text(encoding="utf-8"))

    def test_result_kind_sources_and_cycle_definition(self):
        payload = self.payload
        self.assertEqual(
            payload["result_kind"],
            "SCRATCH_STRICT_COMPLETION_PHASE_SIGN_CYCLE_CLOSURE_OBSTRUCTION_CERTIFICATE__GF2_BOUNDARY_CONDITION",
        )
        self.assertIn("cycle-closure-parity", payload["status"])
        self.assertEqual(set(payload["source_reports"]), {
            "phase_sign_z2_coboundary_certificate",
            "phase_sign_path_cohomology_triviality_certificate",
            "phase_sign_gf2_linear_system_certificate",
        })
        self.assertIn("cycle closure parity", payload["grep_disambiguation"]["searched_terms"])
        definition = payload["cycle_complex_definition"]
        self.assertEqual(definition["field"], "GF(2)")
        self.assertEqual(definition["node_count_dim_C0"], 12)
        self.assertEqual(definition["cycle_edge_count_dim_C1"], 12)
        self.assertEqual(definition["cycle_rank_E_minus_V_plus_components"], 1)
        self.assertEqual(definition["edge_labels"][-1], "11->0")
        self.assertEqual(len(definition["cycle_coboundary_delta_matrix"]), 12)

    def test_rank_kernel_closure_cases_and_summary(self):
        payload = self.payload
        rank = payload["cycle_coboundary_rank_certificate"]
        self.assertEqual(rank["rank"], 11)
        self.assertEqual(rank["nullity"], 1)
        self.assertEqual([row["pivot_col"] for row in rank["pivot_rows"]], list(range(11)))

        kernel = payload["kernel_certificate"]
        self.assertEqual(kernel["kernel_size"], 2)
        self.assertTrue(kernel["kernel_is_exactly_constant_node_cochains"])
        self.assertEqual(kernel["kernel_vectors"], [[0] * 12, [1] * 12])

        cases = {row["case"]: row for row in payload["closure_case_rows"]}
        zero_case = cases["audited_path_plus_forced_zero_closing_edge"]
        self.assertEqual(zero_case["closing_edge_bit"], 0)
        self.assertTrue(zero_case["cycle_parity_even"])
        self.assertTrue(zero_case["is_exact_cycle_coboundary"])
        self.assertEqual(zero_case["exact_node_potential_count"], 2)
        self.assertTrue(zero_case["anchored_reconstruction_matches_audited_nodes"])
        self.assertTrue(zero_case["constant_kernel_pair_present_when_exact"])

        odd_case = cases["odd_closing_edge_perturbation"]
        self.assertEqual(odd_case["closing_edge_bit"], 1)
        self.assertFalse(odd_case["cycle_parity_even"])
        self.assertFalse(odd_case["is_exact_cycle_coboundary"])
        self.assertEqual(odd_case["exact_node_potential_count"], 0)
        self.assertFalse(odd_case["closing_edge_matches_audited_endpoint_xor"])

        comparison_rows = payload["cycle_edge_comparison_rows"]
        self.assertEqual(len(comparison_rows), 12)
        self.assertTrue(all(row["differs_only_on_closing_edge"] for row in comparison_rows))

        summary = payload["cycle_closure_summary"]
        self.assertEqual(summary["rank_delta_cycle"], 11)
        self.assertEqual(summary["nullity_delta_cycle"], 1)
        self.assertEqual(summary["h1_dimension_dim_C1_minus_rank_delta"], 1)
        self.assertEqual(summary["cycle_rank"], 1)
        self.assertEqual(summary["forced_closing_edge_bit_b11_xor_b0"], 0)
        self.assertTrue(summary["zero_closing_cycle_exact"])
        self.assertTrue(summary["zero_closing_anchor_recovers_audited_nodes"])
        self.assertFalse(summary["odd_closing_cycle_exact"])
        self.assertTrue(summary["odd_closing_obstructed_by_cycle_parity"])
        self.assertTrue(summary["path_h1_zero_inherited"])
        self.assertTrue(summary["path_anchor_reconstruction_inherited"])
        self.assertTrue(summary["gf2_path_solution_inherited"])

    def test_proof_blockers_and_hard_limits(self):
        payload = self.payload
        blocker = "\n".join(payload["blocker_context"]["still_open"])
        self.assertIn("strict_phase_frequency_derivation_from_nadsoliton_dynamics", blocker)
        self.assertIn("QW-2191_selector_obstruction", blocker)
        self.assertIn("not a cyclic gate-generation strategy", payload["blocker_context"]["guardrail"])

        proof = payload["proof_certificate"]
        self.assertIn("V=12", proof["cycle_rank_step"])
        self.assertIn("total cycle parity is 0", proof["criterion_step"])
        self.assertIn("b(11) xor b(0)=0", proof["zero_closure_step"])
        self.assertIn("no node potential", proof["odd_closure_step"])
        self.assertIn("does not derive omega/phi", proof["theoretical_limit"])

        hard_limits = "\n".join(payload["hard_limits"])
        self.assertIn("K_strict_gate remains the current live/full", hard_limits)
        self.assertIn("No unqualified identity K_legacy_ont == K_strict_gate", hard_limits)
        self.assertIn("No cyclic L5/L12 gate-generation strategy", hard_limits)
        self.assertIn("No QW-2191 selector discharge", hard_limits)
        self.assertIn("No ToE closure", hard_limits)
        self.assertTrue(REPORT_MD.exists())


if __name__ == "__main__":
    unittest.main()
