import json
import subprocess
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "fundamental_action_reconstruction" / "scratch" / "bridge_strict_completion_theorem_frontier_cut_certificate_probe.py"
REPORT_JSON = ROOT / "fundamental_action_reconstruction" / "scratch" / "bridge_strict_completion_theorem_frontier_cut_certificate_report.json"
REPORT_MD = ROOT / "fundamental_action_reconstruction" / "scratch" / "bridge_strict_completion_theorem_frontier_cut_certificate_report.md"


class TheoremFrontierCutCertificateProbeTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        subprocess.run([sys.executable, str(SCRIPT)], check=True, cwd=ROOT)
        cls.payload = json.loads(REPORT_JSON.read_text(encoding="utf-8"))

    def test_result_kind_sources_and_atoms(self):
        payload = self.payload
        self.assertEqual(
            payload["result_kind"],
            "SCRATCH_STRICT_COMPLETION_THEOREM_FRONTIER_CUT_CERTIFICATE__NO_CLOSURE_THEOREM",
        )
        self.assertIn("theorem-frontier-dag", payload["status"])
        self.assertEqual(set(payload["source_reports"]), {
            "finite_bridge_assembly",
            "symbolic_cancellation",
            "role_transfer_minimal_obligation_lattice",
            "anchor_h1_classification",
            "component_gap_matrix",
            "closure_plan_dependency",
        })
        self.assertIn("theorem readiness", payload["grep_disambiguation"]["searched_terms"])
        self.assertIn("remaining theorem-frontier cut", payload["grep_disambiguation"]["finding"])
        self.assertEqual(len(payload["closed_atoms"]), 4)
        self.assertEqual(len(payload["open_atoms"]), 7)

    def test_dag_targets_summary_and_cross_checks(self):
        payload = self.payload
        summary = payload["theorem_frontier_cut_summary"]
        self.assertEqual(summary["dag_node_count"], 15)
        self.assertTrue(summary["dag_is_acyclic"])
        self.assertEqual(len(summary["topological_order"]), 15)
        self.assertEqual(summary["closed_atom_count"], 4)
        self.assertEqual(summary["open_atom_count"], 7)
        self.assertTrue(summary["all_closed_atoms_certified"])
        self.assertTrue(summary["all_open_atoms_still_missing"])
        self.assertEqual(summary["minimal_open_leaf_cut_for_toe_size"], 7)
        self.assertEqual(set(summary["minimal_open_leaf_cut_for_toe"]), set(payload["open_atoms"]))
        self.assertTrue(summary["component_gap_sources_still_missing"])
        self.assertTrue(summary["closure_plan_keeps_toe_open"])
        self.assertFalse(summary["bridge_theorem_exported"])
        self.assertFalse(summary["role_transfer_theorem_exported"])
        self.assertFalse(summary["selector_closure_exported"])
        self.assertFalse(summary["toe_closure_claimed"])

        targets = {row["target"]: row for row in payload["target_frontier_rows"]}
        self.assertEqual(set(targets), {
            "bridge_theorem_level_closure",
            "role_transfer_theorem_level_closure",
            "selector_qw2191_closure",
            "toe_closure",
        })
        self.assertEqual(targets["bridge_theorem_level_closure"]["open_leaf_cut_size"], 3)
        self.assertEqual(targets["role_transfer_theorem_level_closure"]["open_leaf_cut_size"], 4)
        self.assertEqual(targets["selector_qw2191_closure"]["open_leaf_cut"], ["chi11_selector_source"])
        self.assertTrue(all(payload["cross_checks"].values()))

    def test_proof_and_hard_limits(self):
        payload = self.payload
        proof = payload["proof_certificate"]
        self.assertIn("rg was used", proof["grep_step"])
        self.assertIn("15 nodes", proof["dag_step"])
        self.assertIn("closed certificate atoms", proof["closed_step"])
        self.assertIn("seven missing theorem atoms", proof["open_cut_step"])
        self.assertIn("exports no bridge theorem", proof["scope_step"])

        hard_limits = "\n".join(payload["hard_limits"])
        self.assertIn("No missing theorem atom", hard_limits)
        self.assertIn("No bridge theorem-level closure", hard_limits)
        self.assertIn("No legacy role-transfer theorem", hard_limits)
        self.assertIn("No selector/QW-2191 closure", hard_limits)
        self.assertIn("No ToE closure", hard_limits)
        self.assertTrue(REPORT_MD.exists())


if __name__ == "__main__":
    unittest.main()
