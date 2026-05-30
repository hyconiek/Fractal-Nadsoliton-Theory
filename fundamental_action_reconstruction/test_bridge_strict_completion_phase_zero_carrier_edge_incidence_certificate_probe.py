from __future__ import annotations

import json
import subprocess
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "fundamental_action_reconstruction" / "scratch" / "bridge_strict_completion_phase_zero_carrier_edge_incidence_certificate_probe.py"
REPORT_JSON = ROOT / "fundamental_action_reconstruction" / "scratch" / "bridge_strict_completion_phase_zero_carrier_edge_incidence_certificate_report.json"
REPORT_MD = ROOT / "fundamental_action_reconstruction" / "scratch" / "bridge_strict_completion_phase_zero_carrier_edge_incidence_certificate_report.md"


class StrictCompletionPhaseZeroCarrierEdgeIncidenceCertificateProbeTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        subprocess.run([sys.executable, str(SCRIPT)], check=True, cwd=ROOT)
        cls.payload = json.loads(REPORT_JSON.read_text(encoding="utf-8"))

    def test_result_kind_sources_and_definition(self):
        payload = self.payload
        self.assertEqual(
            payload["result_kind"],
            "SCRATCH_STRICT_COMPLETION_PHASE_ZERO_CARRIER_EDGE_INCIDENCE_CERTIFICATE__GF2_CARRIER_TO_EDGE_PARITY",
        )
        self.assertIn("incidence", payload["status"])
        self.assertEqual(set(payload["source_reports"]), {
            "cell_partition_certificate",
            "cell_sign_certificate",
            "phase_sign_z2_coboundary_certificate",
            "phase_sign_gf2_linear_system_certificate",
        })
        self.assertIn("carrier_edge_incidence_matrix", payload["grep_disambiguation"]["searched_terms"])
        definition = payload["incidence_definition"]
        self.assertEqual(definition["field"], "GF(2)")
        self.assertEqual(definition["carrier_order"], ["legacy_z0", "legacy_z1", "strict_z0", "legacy_z2"])
        self.assertEqual(definition["carrier_multiplicity_vector"], [1, 1, 1, 1])
        self.assertIn("M * carrier_multiplicity_vector", definition["edge_bit_rule"])

    def test_containment_matrix_rank_and_rows(self):
        payload = self.payload
        containment = payload["carrier_rational_containment_rows"]
        self.assertEqual(len(containment), 4)
        self.assertTrue(all(row["strictly_inside_open_edge_by_rational_bounds"] for row in containment))
        self.assertEqual([row["edge"] for row in containment], ["1->2", "5->6", "7->8", "9->10"])

        matrix = payload["carrier_edge_incidence_matrix"]
        self.assertEqual(len(matrix), 11)
        self.assertTrue(all(len(row) == 4 for row in matrix))
        self.assertEqual([row[0] for row in matrix], [0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0])
        self.assertEqual([row[1] for row in matrix], [0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0])
        self.assertEqual([row[2] for row in matrix], [0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0])
        self.assertEqual([row[3] for row in matrix], [0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0])

        columns = payload["carrier_column_rows"]
        self.assertTrue(all(row["is_single_edge_standard_basis_column"] for row in columns))
        self.assertTrue(all(row["matches_declared_carrier_edge"] for row in columns))

        rank = payload["rank_certificate"]
        self.assertEqual(rank["column_rank_over_gf2"], 4)
        self.assertTrue(rank["full_column_rank_on_carrier_subspace"])

    def test_edge_summary_proof_blockers_and_limits(self):
        payload = self.payload
        edge_rows = payload["edge_incidence_rows"]
        self.assertEqual(len(edge_rows), 11)
        self.assertEqual([row["edge"] for row in edge_rows if row["edge_bit_from_incidence"]], ["1->2", "5->6", "7->8", "9->10"])
        self.assertTrue(all(row["edge_bit_equals_odd_carrier_parity"] for row in edge_rows))
        self.assertLessEqual(max(row["incident_carrier_count"] for row in edge_rows), 1)

        summary = payload["carrier_edge_incidence_summary"]
        self.assertEqual(summary["edge_bit_pattern_from_carrier_incidence"], [0, 1, 0, 0, 0, 1, 0, 1, 0, 1, 0])
        self.assertEqual(summary["derived_flip_edges_from_carrier_incidence"], ["1->2", "5->6", "7->8", "9->10"])
        self.assertEqual(summary["column_rank_over_gf2"], 4)
        self.assertTrue(summary["all_carriers_strictly_inside_open_edges_by_rational_bounds"])
        self.assertTrue(summary["all_carrier_columns_are_single_edge_standard_basis_columns"])
        self.assertTrue(summary["all_edges_have_at_most_one_incident_carrier"])
        self.assertTrue(summary["all_edge_bits_equal_odd_carrier_parity"])
        self.assertTrue(summary["matches_expected_carrier_order"])
        self.assertTrue(summary["matches_expected_edge_bits"])
        self.assertTrue(summary["matches_expected_flip_edges"])
        self.assertTrue(summary["matches_cell_partition_flip_edges"])
        self.assertTrue(summary["matches_cell_sign_flip_edges"])
        self.assertTrue(summary["matches_z2_edge_bits"])
        self.assertTrue(summary["matches_gf2_solution_edge_bits"])

        blocker = "\n".join(payload["blocker_context"]["still_open"])
        self.assertIn("strict_phase_frequency_derivation_from_nadsoliton_dynamics", blocker)
        self.assertIn("QW-2191_selector_obstruction", blocker)

        proof = payload["proof_certificate"]
        self.assertIn("positive rational clearance", proof["containment_step"])
        self.assertIn("incidence matrix", proof["incidence_step"])
        self.assertIn("all-one carrier multiplicity vector", proof["parity_step"])
        self.assertIn("column rank 4", proof["rank_step"])
        self.assertIn("does not derive zero carriers", proof["theoretical_limit"])

        hard_limits = "\n".join(payload["hard_limits"])
        self.assertIn("K_strict_gate remains the current live/full", hard_limits)
        self.assertIn("No unqualified identity K_legacy_ont == K_strict_gate", hard_limits)
        self.assertIn("No beta_tors -> chi_11 theorem", hard_limits)
        self.assertIn("No QW-2191 selector discharge", hard_limits)
        self.assertIn("No ToE closure", hard_limits)
        self.assertTrue(REPORT_MD.exists())


if __name__ == "__main__":
    unittest.main()
