from __future__ import annotations

import json
import subprocess
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "fundamental_action_reconstruction" / "scratch" / "bridge_strict_completion_phase_zero_carrier_prefix_node_matrix_certificate_probe.py"
REPORT_JSON = ROOT / "fundamental_action_reconstruction" / "scratch" / "bridge_strict_completion_phase_zero_carrier_prefix_node_matrix_certificate_report.json"
REPORT_MD = ROOT / "fundamental_action_reconstruction" / "scratch" / "bridge_strict_completion_phase_zero_carrier_prefix_node_matrix_certificate_report.md"


class StrictCompletionPhaseZeroCarrierPrefixNodeMatrixCertificateProbeTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        subprocess.run([sys.executable, str(SCRIPT)], check=True, cwd=ROOT)
        cls.payload = json.loads(REPORT_JSON.read_text(encoding="utf-8"))

    def test_result_kind_sources_and_definition(self):
        payload = self.payload
        self.assertEqual(
            payload["result_kind"],
            "SCRATCH_STRICT_COMPLETION_PHASE_ZERO_CARRIER_PREFIX_NODE_MATRIX_CERTIFICATE__GF2_CARRIER_TO_NODE_PARITY",
        )
        self.assertIn("carrier-prefix", payload["status"])
        self.assertEqual(set(payload["source_reports"]), {
            "carrier_edge_incidence_certificate",
            "cell_sign_certificate",
            "phase_sign_z2_coboundary_certificate",
            "phase_sign_gf2_linear_system_certificate",
        })
        self.assertIn("carrier-prefix", payload["grep_disambiguation"]["searched_terms"])
        definition = payload["matrix_definition"]
        self.assertEqual(definition["field"], "GF(2)")
        self.assertEqual(definition["carrier_order"], ["legacy_z0", "legacy_z1", "strict_z0", "legacy_z2"])
        self.assertEqual(definition["carrier_multiplicity_vector"], [1, 1, 1, 1])
        self.assertEqual(len(definition["carrier_prefix_full_node_matrix_C"]), 12)
        self.assertIn("C_tail * carrier_multiplicity_vector", definition["node_bit_rule"])

    def test_node_rows_edge_differences_and_rank(self):
        payload = self.payload
        rows = payload["node_prefix_rows"]
        self.assertEqual(len(rows), 12)
        self.assertEqual(rows[0]["carriers_left_by_matrix"], [])
        self.assertEqual(rows[2]["carriers_left_by_matrix"], ["legacy_z0"])
        self.assertEqual(rows[8]["carriers_left_by_matrix"], ["legacy_z0", "legacy_z1", "strict_z0"])
        self.assertEqual(rows[10]["carriers_left_by_matrix"], ["legacy_z0", "legacy_z1", "strict_z0", "legacy_z2"])
        self.assertTrue(all(row["matches_cell_sign_left_carriers"] for row in rows))
        self.assertTrue(all(row["matches_cell_sign_node_bit"] for row in rows))
        self.assertEqual([row["node_bit_from_matrix"] for row in rows], [0, 0, 1, 1, 1, 1, 0, 0, 1, 1, 0, 0])

        edge_rows = payload["edge_difference_rows"]
        self.assertEqual(len(edge_rows), 11)
        self.assertTrue(all(row["difference_recovers_edge_incidence"] for row in edge_rows))
        self.assertEqual(edge_rows[1]["node_prefix_difference_row"], [1, 0, 0, 0])
        self.assertEqual(edge_rows[7]["node_prefix_difference_row"], [0, 0, 1, 0])

        rank = payload["rank_certificate"]
        self.assertEqual(rank["carrier_prefix_rank_over_gf2"], 4)
        self.assertTrue(rank["full_column_rank_on_carrier_prefix_subspace"])

    def test_summary_proof_blockers_and_limits(self):
        payload = self.payload
        summary = payload["carrier_prefix_node_matrix_summary"]
        self.assertEqual(summary["carrier_order"], ["legacy_z0", "legacy_z1", "strict_z0", "legacy_z2"])
        self.assertEqual(summary["node_bit_pattern_from_carrier_prefix_matrix"], [0, 0, 1, 1, 1, 1, 0, 0, 1, 1, 0, 0])
        self.assertEqual(summary["phase_sign_pattern_from_carrier_prefix_matrix"], [1, 1, -1, -1, -1, -1, 1, 1, -1, -1, 1, 1])
        self.assertEqual(summary["carrier_prefix_rank_over_gf2"], 4)
        self.assertTrue(summary["all_node_rows_match_cell_sign_left_carriers"])
        self.assertTrue(summary["all_node_bits_match_cell_sign"])
        self.assertTrue(summary["all_edge_differences_recover_carrier_edge_incidence"])
        self.assertTrue(summary["matches_expected_carrier_order"])
        self.assertTrue(summary["matches_expected_node_bits"])
        self.assertTrue(summary["matches_expected_sign_pattern"])
        self.assertTrue(summary["matches_z2_node_bits"])
        self.assertTrue(summary["matches_gf2_prefix_solution"])
        self.assertTrue(summary["inherits_carrier_edge_incidence_rank"])

        blocker = "\n".join(payload["blocker_context"]["still_open"])
        self.assertIn("strict_phase_frequency_derivation_from_nadsoliton_dynamics", blocker)
        self.assertIn("QW-2191_selector_obstruction", blocker)

        proof = payload["proof_certificate"]
        self.assertIn("C=L*M", proof["composition_step"])
        self.assertIn("audited node-bit pattern", proof["node_step"])
        self.assertIn("carriers-left-of-node", proof["cell_sign_step"])
        self.assertIn("Adjacent node-prefix row differences", proof["edge_difference_step"])
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
