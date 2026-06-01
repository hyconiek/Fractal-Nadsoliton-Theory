import json
import subprocess
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "fundamental_action_reconstruction" / "scratch" / "bridge_strict_completion_theorem_frontier_truth_table_certificate_probe.py"
REPORT_JSON = ROOT / "fundamental_action_reconstruction" / "scratch" / "bridge_strict_completion_theorem_frontier_truth_table_certificate_report.json"
REPORT_MD = ROOT / "fundamental_action_reconstruction" / "scratch" / "bridge_strict_completion_theorem_frontier_truth_table_certificate_report.md"


class TheoremFrontierTruthTableCertificateProbeTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        subprocess.run([sys.executable, str(SCRIPT)], check=True, cwd=ROOT)
        cls.payload = json.loads(REPORT_JSON.read_text(encoding="utf-8"))

    def test_result_kind_sources_and_atoms(self):
        payload = self.payload
        self.assertEqual(
            payload["result_kind"],
            "SCRATCH_STRICT_COMPLETION_THEOREM_FRONTIER_TRUTH_TABLE_CERTIFICATE__NO_CLOSURE_THEOREM",
        )
        self.assertIn("all-2pow7-frontier-assignments", payload["status"])
        self.assertEqual(set(payload["source_reports"]), {
            "theorem_frontier_cut",
            "role_transfer_minimal_obligation_lattice",
            "component_gap_matrix",
            "anchor_h1_classification",
        })
        self.assertEqual(len(payload["open_atoms"]), 7)
        self.assertEqual(payload["target_definitions"]["bridge_theorem_level_closure"], [
            "strict_dynamical_source_for_A_P_D",
            "strict_phase_frequency_source",
            "strict_damping_beta_eta_source",
        ])
        self.assertEqual(payload["target_definitions"]["selector_qw2191_closure"], ["chi11_selector_source"])
        self.assertIn("truth table", payload["grep_disambiguation"]["searched_terms"])
        self.assertIn("finite Boolean audit", payload["grep_disambiguation"]["finding"])

    def test_truth_table_summary_and_targets(self):
        payload = self.payload
        summary = payload["theorem_frontier_truth_table_summary"]
        self.assertEqual(summary["open_atom_count"], 7)
        self.assertEqual(summary["truth_assignment_count"], 128)
        self.assertTrue(summary["current_assignment_all_false"])
        self.assertTrue(summary["current_targets_all_false"])
        self.assertEqual(summary["bridge_satisfying_assignment_count"], 16)
        self.assertEqual(summary["role_satisfying_assignment_count"], 8)
        self.assertEqual(summary["selector_satisfying_assignment_count"], 64)
        self.assertEqual(summary["toe_satisfying_assignment_count"], 1)
        self.assertEqual(summary["bridge_minimal_set_size"], 3)
        self.assertEqual(summary["role_minimal_set_size"], 4)
        self.assertEqual(summary["selector_minimal_set_size"], 1)
        self.assertEqual(summary["toe_minimal_set_size"], 7)
        self.assertTrue(summary["toe_minimal_set_equals_frontier_cut"])
        self.assertTrue(summary["role_lattice_minimal_set_inherited"])
        self.assertTrue(summary["component_gap_sources_still_missing"])
        self.assertTrue(summary["anchor_selector_source_still_open"])
        self.assertFalse(summary["bridge_theorem_exported"])
        self.assertFalse(summary["role_transfer_theorem_exported"])
        self.assertFalse(summary["selector_closure_exported"])
        self.assertFalse(summary["toe_closure_claimed"])

        self.assertEqual(len(payload["truth_table_rows"]), 128)
        self.assertEqual(payload["truth_table_rows"][0]["true_atoms"], [])
        self.assertFalse(payload["truth_table_rows"][0]["toe_closure"])
        self.assertEqual(payload["truth_table_rows"][-1]["true_atom_count"], 7)
        self.assertTrue(payload["truth_table_rows"][-1]["toe_closure"])
        self.assertEqual(payload["target_satisfying_assignment_counts"], {
            "bridge_theorem_level_closure": 16,
            "role_transfer_theorem_level_closure": 8,
            "selector_qw2191_closure": 64,
            "toe_closure": 1,
        })
        self.assertEqual(payload["target_minimal_true_atom_sets"]["toe_closure"], [payload["open_atoms"]])
        self.assertTrue(all(payload["cross_checks"].values()))

    def test_proof_and_hard_limits(self):
        payload = self.payload
        proof = payload["proof_certificate"]
        self.assertIn("rg was used", proof["grep_step"])
        self.assertIn("All 2^7=128 truth assignments", proof["enumeration_step"])
        self.assertIn("Bridge closure requires three", proof["target_step"])
        self.assertIn("all seven atoms false", proof["current_step"])
        self.assertIn("bridge=3, role=4, selector=1, and ToE=7", proof["minimal_step"])
        self.assertIn("exports no theorem atom", proof["scope_step"])

        hard_limits = "\n".join(payload["hard_limits"])
        self.assertIn("No truth-table assignment", hard_limits)
        self.assertIn("No missing theorem atom", hard_limits)
        self.assertIn("No bridge theorem", hard_limits)
        self.assertIn("No QW-2191 selector discharge", hard_limits)
        self.assertTrue(REPORT_MD.exists())


if __name__ == "__main__":
    unittest.main()
