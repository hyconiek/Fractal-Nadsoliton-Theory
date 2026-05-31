from __future__ import annotations

import json
import subprocess
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "fundamental_action_reconstruction" / "scratch" / "bridge_strict_completion_phase_sign_node_support_interval_boundary_certificate_probe.py"
REPORT_JSON = ROOT / "fundamental_action_reconstruction" / "scratch" / "bridge_strict_completion_phase_sign_node_support_interval_boundary_certificate_report.json"
REPORT_MD = ROOT / "fundamental_action_reconstruction" / "scratch" / "bridge_strict_completion_phase_sign_node_support_interval_boundary_certificate_report.md"


class StrictCompletionPhaseSignNodeSupportIntervalBoundaryCertificateProbeTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        subprocess.run([sys.executable, str(SCRIPT)], check=True, cwd=ROOT)
        cls.payload = json.loads(REPORT_JSON.read_text(encoding="utf-8"))

    def test_result_kind_sources_and_interval_definition(self):
        payload = self.payload
        self.assertEqual(
            payload["result_kind"],
            "SCRATCH_STRICT_COMPLETION_PHASE_SIGN_NODE_SUPPORT_INTERVAL_BOUNDARY_CERTIFICATE__GF2_COMPONENT_BOUNDARY",
        )
        self.assertIn("node-support-maximal-interval-boundary", payload["status"])
        self.assertEqual(set(payload["source_reports"]), {
            "phase_sign_z2_coboundary_certificate",
            "phase_sign_edge_support_minimality_certificate",
            "phase_sign_reduced_coboundary_inverse_certificate",
            "phase_sign_cycle_closure_obstruction_certificate",
        })
        self.assertIn("component boundary", payload["grep_disambiguation"]["searched_terms"])
        definition = payload["node_support_interval_definition"]
        self.assertEqual(definition["field"], "GF(2)")
        self.assertEqual(definition["one_node_support"], [2, 3, 4, 5, 8, 9])
        self.assertEqual(definition["maximal_one_intervals"], [{"start": 2, "end": 5}, {"start": 8, "end": 9}])
        self.assertEqual(definition["interval_multiplicity_vector"], [1, 1])
        self.assertEqual(len(definition["interval_boundary_matrix_edges_by_intervals"]), 11)

    def test_interval_rows_rank_boundary_formula_and_reconstruction(self):
        payload = self.payload
        rows = payload["interval_rows"]
        self.assertEqual(len(rows), 2)
        self.assertEqual(rows[0]["boundary_edges"], ["1->2", "5->6"])
        self.assertEqual(rows[1]["boundary_edges"], ["7->8", "9->10"])
        self.assertTrue(all(not row["touches_left_endpoint"] and not row["touches_right_endpoint"] for row in rows))
        self.assertEqual([row["boundary_weight"] for row in rows], [2, 2])

        rank = payload["rank_certificate"]
        self.assertEqual(rank["interval_boundary_column_rank_over_gf2"], 2)
        self.assertEqual(rank["interval_boundary_column_nullity_over_gf2"], 0)
        self.assertTrue(rank["full_column_rank_on_interval_boundary_subspace"])

        formula = payload["component_boundary_formula"]
        self.assertEqual(formula["component_count"], 2)
        self.assertEqual(formula["endpoint_touch_count"], 0)
        self.assertEqual(formula["predicted_boundary_weight_2c_minus_endpoint_touches"], 4)
        self.assertEqual(formula["actual_edge_bit_hamming_weight"], 4)
        self.assertTrue(formula["formula_matches"])

        reconstruction = payload["reconstruction_certificate"]
        self.assertEqual(reconstruction["node_bits_from_interval_support_union"], [0, 0, 1, 1, 1, 1, 0, 0, 1, 1, 0, 0])
        self.assertEqual(reconstruction["edge_bits_from_interval_boundaries"], [0, 1, 0, 0, 0, 1, 0, 1, 0, 1, 0])
        self.assertEqual(reconstruction["flip_edges_from_interval_boundaries"], ["1->2", "5->6", "7->8", "9->10"])
        self.assertTrue(all(row["matches_audited_edge_bit"] for row in payload["edge_boundary_rows"]))

    def test_summary_proof_blockers_and_hard_limits(self):
        payload = self.payload
        summary = payload["node_support_interval_boundary_summary"]
        self.assertEqual(summary["maximal_interval_count"], 2)
        self.assertTrue(summary["matches_expected_intervals"])
        self.assertTrue(summary["interval_boundary_full_column_rank"])
        self.assertTrue(summary["node_bits_recovered_from_interval_union"])
        self.assertTrue(summary["edge_bits_recovered_from_interval_boundaries"])
        self.assertTrue(summary["boundary_weight_formula_matches"])
        self.assertTrue(summary["matches_edge_support_minimality_flip_edges"])
        self.assertTrue(summary["matches_reduced_coboundary_edge_bits"])
        self.assertTrue(summary["contrasts_cycle_odd_closing_obstruction"])

        blocker = "\n".join(payload["blocker_context"]["still_open"])
        self.assertIn("strict_phase_frequency_derivation_from_nadsoliton_dynamics", blocker)
        self.assertIn("QW-2191_selector_obstruction", blocker)

        proof = payload["proof_certificate"]
        self.assertIn("[2,5] U [8,9]", proof["support_step"])
        self.assertIn("edge bits [0,1,0,0,0,1,0,1,0,1,0]", proof["boundary_step"])
        self.assertIn("rank 2", proof["rank_step"])
        self.assertIn("2*component_count = 4", proof["component_count_step"])
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
