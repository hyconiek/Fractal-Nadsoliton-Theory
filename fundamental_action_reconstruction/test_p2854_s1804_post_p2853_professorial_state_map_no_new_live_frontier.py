import json
import subprocess
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
P2853_SCRIPT = ROOT / "p2853_s1803_phase_frequency_bridge_source_audit.py"
SCRIPT = ROOT / "p2854_s1804_post_p2853_professorial_state_map_no_new_live_frontier.py"
JSON_PATH = ROOT / "generated" / "p2854_s1804_post_p2853_professorial_state_map_no_new_live_frontier.json"
MD_PATH = ROOT / "generated" / "p2854_s1804_post_p2853_professorial_state_map_no_new_live_frontier.md"


class P2854PostP2853ProfessorialStateMapTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        subprocess.run([sys.executable, str(P2853_SCRIPT)], check=True, cwd=ROOT)
        subprocess.run([sys.executable, str(SCRIPT)], check=True, cwd=ROOT)
        cls.payload = json.loads(JSON_PATH.read_text(encoding="utf-8"))
        cls.state = cls.payload["post_p2853_state_map"]

    def test_status_and_inputs(self):
        self.assertEqual(self.payload["status"], "P2854_POST_P2853_PROFESSORIAL_STATE_MAP_NO_NEW_LIVE_FRONTIER")
        statuses = self.state["input_statuses_rechecked"]
        self.assertEqual(statuses["P2848"], "P2848_COUPLING_COEFFICIENT_UNIT_SOURCE_LAW_AUDIT_NO_CLOSURE")
        self.assertEqual(statuses["P2849"], "P2849_DAMPING_COMPRESSION_KERNEL_BRIDGE_ATOM_AUDIT_NO_CLOSURE")
        self.assertEqual(statuses["P2852"], "P2852_KERNEL_BRIDGE_OBLIGATION_RECONCILIATION_MATRIX_NO_CLOSURE")
        self.assertEqual(statuses["P2853"], "P2853_PHASE_FREQUENCY_BRIDGE_SOURCE_AUDIT_NO_CLOSURE")

    def test_no_live_frontier_counts(self):
        self.assertEqual(self.state["closure_row_count"], 0)
        self.assertEqual(self.state["admissible_without_new_premise_count"], 0)
        self.assertTrue(self.payload["acceptance_matrix"]["accepted_as_no_new_live_frontier_certificate"])
        self.assertFalse(self.payload["acceptance_matrix"]["exports_new_source_premise"])

    def test_frontier_rows_cover_expected_lanes(self):
        rows = {row["lane"]: row for row in self.state["frontier_rows"]}
        for lane in (
            "finite_density_to_unit_bearing_ltotal",
            "damping_compression_eta_beta",
            "eml_syntax_basis",
            "amplitude_alpha_geo_passage",
            "phase_frequency_omega_phi_source",
            "full_kernel_completion_bridge",
            "selector_topological_qw2191",
            "legacy_role_transfer",
            "ltotal_eom_hamiltonian_toe",
        ):
            self.assertIn(lane, rows)
            self.assertFalse(rows[lane]["closure_exported_now"])
            self.assertFalse(rows[lane]["admissible_without_new_premise"])
        self.assertIn("omega=743/4000", rows["phase_frequency_omega_phi_source"]["missing_new_premise"])
        self.assertIn("target-independent positive beta source", rows["damping_compression_eta_beta"]["missing_new_premise"])

    def test_professorial_path_and_source_tests(self):
        tests = {row["candidate_source"]: row for row in self.state["proof_grade_source_acceptance_tests"]}
        self.assertIn("strict_phase_frequency_source_law", tests)
        self.assertIn("strict_eta_beta_source_law", tests)
        self.assertIn("role_safe_alpha_geo_source_law", tests)
        self.assertIn("omega=743/4000", tests["strict_phase_frequency_source_law"]["first_test"])
        gates = [row["gate"] for row in self.state["professorial_closure_path"]]
        self.assertEqual(gates, ["sourcehood", "component_bridge", "completion_map", "role_transfer", "action_eom_hamiltonian"])

    def test_no_false_closure_and_documents(self):
        flags = self.payload["decision"]["negative_export_flags"]
        self.assertFalse(flags["new_live_frontier_without_new_premise_exported"])
        self.assertFalse(flags["strict_phase_frequency_source_law_exported"])
        self.assertFalse(flags["eta_beta_source_law_exported"])
        self.assertFalse(flags["full_kernel_bridge_exported"])
        self.assertFalse(flags["ltotal_exported"])
        self.assertFalse(flags["toe_closure_exported"])
        self.assertIn("P2854/S1804", MD_PATH.read_text(encoding="utf-8"))
        self.assertIn("P2854/S1804", (ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md").read_text(encoding="utf-8"))
        self.assertIn("P2854/S1804", (ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md").read_text(encoding="utf-8"))
        self.assertIn("P2854", (REPO / "AGENTS.md").read_text(encoding="utf-8"))


if __name__ == "__main__":
    unittest.main()
