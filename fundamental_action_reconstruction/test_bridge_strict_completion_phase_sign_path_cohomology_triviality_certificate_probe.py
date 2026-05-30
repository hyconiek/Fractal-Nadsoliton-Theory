from __future__ import annotations

import json
import subprocess
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "fundamental_action_reconstruction" / "scratch" / "bridge_strict_completion_phase_sign_path_cohomology_triviality_certificate_probe.py"
REPORT_JSON = ROOT / "fundamental_action_reconstruction" / "scratch" / "bridge_strict_completion_phase_sign_path_cohomology_triviality_certificate_report.json"
REPORT_MD = ROOT / "fundamental_action_reconstruction" / "scratch" / "bridge_strict_completion_phase_sign_path_cohomology_triviality_certificate_report.md"


class StrictCompletionPhaseSignPathCohomologyTrivialityCertificateProbeTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        subprocess.run([sys.executable, str(SCRIPT)], check=True, cwd=ROOT)
        cls.payload = json.loads(REPORT_JSON.read_text(encoding="utf-8"))

    def test_result_kind_sources_and_complex_definition(self):
        payload = self.payload
        self.assertEqual(
            payload["result_kind"],
            "SCRATCH_STRICT_COMPLETION_PHASE_SIGN_PATH_COHOMOLOGY_TRIVIALITY_CERTIFICATE__GF2_NO_CYCLE_AMBIGUITY",
        )
        self.assertIn("h1-zero", payload["status"])
        self.assertEqual(set(payload["source_reports"]), {
            "phase_sign_z2_coboundary_certificate",
            "phase_sign_gf2_linear_system_certificate",
            "phase_sign_edge_support_minimality_certificate",
            "phase_zero_gf2_commutative_diagram_certificate",
        })
        self.assertIn("path cohomology", payload["grep_disambiguation"]["searched_terms"])
        definition = payload["cochain_complex_definition"]
        self.assertEqual(definition["field"], "GF(2)")
        self.assertEqual(definition["node_count_dim_C0"], 12)
        self.assertEqual(definition["edge_count_dim_C1"], 11)
        self.assertEqual(definition["cycle_rank_E_minus_V_plus_components"], 0)
        self.assertEqual(len(definition["coboundary_delta_matrix"]), 11)

    def test_rank_kernel_exactness_rows_and_summary(self):
        payload = self.payload
        rank = payload["coboundary_rank_certificate"]
        self.assertEqual(rank["rank"], 11)
        self.assertEqual(rank["nullity"], 1)
        self.assertEqual([row["pivot_col"] for row in rank["pivot_rows"]], list(range(11)))

        kernel = payload["kernel_certificate"]
        self.assertEqual(kernel["kernel_size"], 2)
        self.assertTrue(kernel["kernel_is_exactly_constant_node_cochains"])
        self.assertEqual(kernel["kernel_vectors"], [[0] * 12, [1] * 12])

        exact_rows = payload["edge_exactness_rows"]
        self.assertEqual(len(exact_rows), 11)
        self.assertTrue(all(row["matches_audited_edge_bit"] for row in exact_rows))
        self.assertEqual([row["edge"] for row in exact_rows if row["is_flip_edge"]], ["1->2", "5->6", "7->8", "9->10"])

        summary = payload["path_cohomology_summary"]
        self.assertEqual(summary["rank_delta"], 11)
        self.assertEqual(summary["nullity_delta"], 1)
        self.assertEqual(summary["h1_dimension_dim_C1_minus_rank_delta"], 0)
        self.assertEqual(summary["cycle_rank"], 0)
        self.assertTrue(summary["every_edge_cochain_exact_on_path"])
        self.assertTrue(summary["anchor_kills_constant_kernel"])
        self.assertEqual(summary["anchored_reconstructed_node_bits"], [0, 0, 1, 1, 1, 1, 0, 0, 1, 1, 0, 0])
        self.assertTrue(summary["anchored_reconstruction_matches_z2_node_bits"])
        self.assertTrue(summary["delta_node_bits_match_edge_bits"])
        self.assertTrue(summary["matches_gf2_linear_system_solution"])
        self.assertTrue(summary["inherits_edge_support_unique_minimality"])
        self.assertTrue(summary["inherits_commutative_diagram"])

    def test_proof_blockers_and_hard_limits(self):
        payload = self.payload
        blocker = "\n".join(payload["blocker_context"]["still_open"])
        self.assertIn("strict_phase_frequency_derivation_from_nadsoliton_dynamics", blocker)
        self.assertIn("QW-2191_selector_obstruction", blocker)

        proof = payload["proof_certificate"]
        self.assertIn("rank 11", proof["rank_step"])
        self.assertIn("constant node cochains", proof["kernel_step"])
        self.assertIn("dimension 0", proof["h1_step"])
        self.assertIn("b(0)=0", proof["anchor_step"])
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
