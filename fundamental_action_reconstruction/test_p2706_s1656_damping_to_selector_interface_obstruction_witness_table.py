from __future__ import annotations

import json
import subprocess
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
SCRIPT = ROOT / "p2706_s1656_damping_to_selector_interface_obstruction_witness_table.py"
OUT = ROOT / "generated" / "p2706_s1656_damping_to_selector_interface_obstruction_witness_table.json"
MD = ROOT / "generated" / "p2706_s1656_damping_to_selector_interface_obstruction_witness_table.md"


class P2706DampingToSelectorInterfaceObstructionTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        subprocess.run([sys.executable, str(SCRIPT)], cwd=ROOT, check=True)
        cls.payload = json.loads(OUT.read_text(encoding="utf-8"))

    def test_finite_orientation_blindness_counts(self):
        count_scan = self.payload["count_scan"]
        self.assertEqual(count_scan["support_count"], 792)
        self.assertEqual(count_scan["worst_directed_count_delta_between_u_and_minus_u"], 0)
        self.assertEqual(count_scan["nonzero_delta_witness_examples"], [])

    def test_weighted_p2377_p2378_scores_are_orientation_blind(self):
        weighted = self.payload["weighted_scan"]
        self.assertEqual(weighted["total_weighted_cases"], 792 * 3 * 3 * 5)
        self.assertEqual(weighted["max_abs_score_delta_plus_direction_vs_minus_direction"], 0.0)
        self.assertEqual(weighted["nonzero_score_delta_examples"], [])

    def test_decision_preserves_no_unlock_boundary(self):
        self.assertEqual(self.payload["status"], "P2706_DAMPING_TO_SELECTOR_INTERFACE_OBSTRUCTION_NO_UNLOCK")
        self.assertTrue(all(row["passes"] for row in self.payload["obstruction_matrix"]))
        decision = self.payload["decision"]
        self.assertFalse(decision["damping_transport_exports_directed_selector"])
        self.assertTrue(all(value is False for value in decision["negative_export_flags"].values()))
        self.assertIn("no-new-live-frontier", decision["next_honest_step"])

    def test_outputs_and_docs_written(self):
        self.assertIn("P2706/S1656", MD.read_text(encoding="utf-8"))
        self.assertIn("P2706/S1656", (ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md").read_text(encoding="utf-8"))
        self.assertIn("P2706/S1656", (ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md").read_text(encoding="utf-8"))
        self.assertIn("P2706", (REPO / "AGENTS.md").read_text(encoding="utf-8"))


if __name__ == "__main__":
    unittest.main()
