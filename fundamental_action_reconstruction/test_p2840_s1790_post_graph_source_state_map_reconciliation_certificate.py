import json
import subprocess
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
SCRIPT = ROOT / "p2840_s1790_post_graph_source_state_map_reconciliation_certificate.py"
JSON_PATH = ROOT / "generated" / "p2840_s1790_post_graph_source_state_map_reconciliation_certificate.json"
MD_PATH = ROOT / "generated" / "p2840_s1790_post_graph_source_state_map_reconciliation_certificate.md"


class P2840PostGraphSourceStateMapReconciliationCertificateTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        subprocess.run([sys.executable, str(SCRIPT)], check=True, cwd=ROOT)
        cls.payload = json.loads(JSON_PATH.read_text(encoding="utf-8"))

    def test_status(self):
        self.assertEqual(
            self.payload["status"],
            "P2840_POST_GRAPH_SOURCE_STATE_MAP_NO_NEW_LIVE_FRONTIER_CERTIFICATE_NO_CLOSURE",
        )

    def test_matrix_counts_and_certificate(self):
        audit = self.payload["post_graph_source_state_map_reconciliation_certificate"]
        self.assertEqual(audit["graph_source_lane_row_count"], 6)
        self.assertEqual(audit["adjacent_replay_gate_row_count"], 5)
        self.assertEqual(audit["failed_evidence_rows"], [])
        self.assertEqual(audit["live_graph_source_frontier_rows"], [])
        self.assertEqual(audit["new_adjacent_object_rows"], [])
        self.assertTrue(audit["no_new_live_frontier_certificate"])

    def test_graph_source_lane_rows(self):
        audit = self.payload["post_graph_source_state_map_reconciliation_certificate"]
        rows = {row["lane"]: row for row in audit["graph_source_lane_rows"]}
        self.assertEqual(set(rows), {
            "finite_combined_separator",
            "target_independent_units_normalization",
            "typed_domain_codomain_map",
            "formal_variational_derivative",
            "evaluation_pullback_localization",
            "unit_bearing_coupling_source_law",
        })
        self.assertEqual(rows["finite_combined_separator"]["status"], "positive_finite_witness_closed")
        self.assertTrue(all(row["evidence_flag"] for row in rows.values()))
        self.assertFalse(any(row["live_frontier_unlocked"] for row in rows.values()))

    def test_adjacent_replay_gates(self):
        audit = self.payload["post_graph_source_state_map_reconciliation_certificate"]
        rows = {row["lane"]: row for row in audit["adjacent_replay_gate_rows"]}
        self.assertEqual(set(rows), {
            "selector_QW2191_replay",
            "legacy_strict_bridge_role_transfer",
            "L_total_or_ToE_promotion",
            "direct_route_residual_replay",
            "lagrangian_EOM_reverse_closure_replay",
        })
        self.assertFalse(any(row["new_object_supplied"] for row in rows.values()))
        self.assertEqual(rows["L_total_or_ToE_promotion"]["status"], "blocked")

    def test_acceptance_boundaries(self):
        matrix = self.payload["acceptance_matrix"]
        facts = matrix["facts"]
        self.assertTrue(facts["all_graph_source_rows_accounted"])
        self.assertTrue(facts["all_adjacent_replay_rows_accounted"])
        self.assertTrue(facts["all_input_evidence_flags_pass"])
        self.assertTrue(facts["no_live_graph_source_frontier_unlocked"])
        self.assertTrue(facts["no_new_adjacent_object_supplied"])
        self.assertFalse(facts["selector_bridge_or_role_transfer_imported"])
        self.assertTrue(matrix["accepted_as_post_graph_source_state_map_reconciliation"])
        self.assertTrue(matrix["accepted_as_no_new_live_frontier_certificate"])
        self.assertFalse(matrix["accepted_as_closure_or_promotion"])

    def test_documents_updated(self):
        self.assertIn("P2840/S1790", MD_PATH.read_text(encoding="utf-8"))
        self.assertIn("P2840/S1790", (ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md").read_text(encoding="utf-8"))
        self.assertIn("P2840/S1790", (ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md").read_text(encoding="utf-8"))
        self.assertIn("P2840", (REPO / "AGENTS.md").read_text(encoding="utf-8"))


if __name__ == "__main__":
    unittest.main()
