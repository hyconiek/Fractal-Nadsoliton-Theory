import json
import subprocess
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
SCRIPT = ROOT / "p2841_s1791_post_p2840_broad_state_map_intake_gate.py"
JSON_PATH = ROOT / "generated" / "p2841_s1791_post_p2840_broad_state_map_intake_gate.json"
MD_PATH = ROOT / "generated" / "p2841_s1791_post_p2840_broad_state_map_intake_gate.md"


class P2841PostP2840BroadStateMapIntakeGateTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        subprocess.run([sys.executable, str(SCRIPT)], check=True, cwd=ROOT)
        cls.payload = json.loads(JSON_PATH.read_text(encoding="utf-8"))

    def test_status(self):
        self.assertEqual(
            self.payload["status"],
            "P2841_POST_P2840_BROAD_STATE_MAP_INTAKE_GATE_NO_NEW_LIVE_FRONTIER_NO_CLOSURE",
        )

    def test_state_map_rows(self):
        audit = self.payload["post_p2840_broad_state_map_intake_gate"]
        rows = {row["lane"]: row for row in audit["state_map_lane_rows"]}
        self.assertEqual(audit["state_map_lane_count"], 8)
        self.assertEqual(set(rows), {
            "graph_source_P2835_P2840",
            "selector_QW2191_orientation",
            "legacy_strict_bridge_role_transfer",
            "direct_route_residual_cancellation",
            "lagrangian_EOM_reverse_closure",
            "lower_boundary_tau_pair_recursion",
            "P2680_non_selector_source_atoms",
            "L_total_ToE_promotion",
        })
        self.assertTrue(all(row["evidence_flag"] for row in rows.values()))
        self.assertEqual(rows["graph_source_P2835_P2840"]["current_status"], "no_new_live_frontier_certificate")

    def test_no_new_object_fields(self):
        audit = self.payload["post_p2840_broad_state_map_intake_gate"]
        self.assertEqual(audit["failed_evidence_rows"], [])
        self.assertEqual(audit["new_object_unlocked_rows"], [])
        self.assertTrue(audit["no_new_live_frontier_certificate"])
        for row in audit["state_map_lane_rows"]:
            self.assertFalse(row["new_strict_typed_object"])
            self.assertFalse(row["new_source_theorem"])
            self.assertFalse(row["new_blocker_cut"])
            self.assertFalse(row["new_provider_class"])
            self.assertFalse(row["new_coupling_or_localization_theorem"])

    def test_acceptance_boundaries(self):
        matrix = self.payload["acceptance_matrix"]
        facts = matrix["facts"]
        self.assertTrue(facts["all_state_map_lanes_accounted"])
        self.assertTrue(facts["all_evidence_flags_pass"])
        self.assertTrue(facts["no_new_typed_object_or_provider_supplied"])
        self.assertTrue(facts["no_new_live_frontier_certificate"])
        self.assertFalse(facts["selector_bridge_or_role_transfer_imported"])
        self.assertTrue(matrix["accepted_as_broad_state_map_intake_gate"])
        self.assertTrue(matrix["accepted_as_no_new_live_frontier_certificate"])
        self.assertFalse(matrix["accepted_as_closure_or_promotion"])

    def test_documents_updated(self):
        self.assertIn("P2841/S1791", MD_PATH.read_text(encoding="utf-8"))
        self.assertIn("P2841/S1791", (ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md").read_text(encoding="utf-8"))
        self.assertIn("P2841/S1791", (ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md").read_text(encoding="utf-8"))
        self.assertIn("P2841", (REPO / "AGENTS.md").read_text(encoding="utf-8"))


if __name__ == "__main__":
    unittest.main()
