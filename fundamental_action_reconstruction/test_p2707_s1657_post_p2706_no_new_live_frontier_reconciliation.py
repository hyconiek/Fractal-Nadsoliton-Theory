from __future__ import annotations

import json
import subprocess
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
SCRIPT = ROOT / "p2707_s1657_post_p2706_no_new_live_frontier_reconciliation.py"
OUT = ROOT / "generated" / "p2707_s1657_post_p2706_no_new_live_frontier_reconciliation.json"
MD = ROOT / "generated" / "p2707_s1657_post_p2706_no_new_live_frontier_reconciliation.md"


class P2707PostP2706NoNewLiveFrontierTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        subprocess.run([sys.executable, str(SCRIPT)], cwd=ROOT, check=True)
        cls.payload = json.loads(OUT.read_text(encoding="utf-8"))

    def test_lane_matrix_has_no_live_frontier(self):
        rows = self.payload["lane_matrix"]
        self.assertGreaterEqual(len(rows), 7)
        self.assertTrue(all(row["passes_no_live_frontier_check"] for row in rows))
        self.assertEqual(self.payload["decision"]["live_lanes_now"], [])

    def test_generated_artifact_inventory_has_no_forbidden_true_flags(self):
        inv = self.payload["post_p2697_inventory"]
        self.assertGreater(inv["p27xx_json_artifacts_scanned"], 0)
        self.assertEqual(inv["artifacts_with_live_frontier_true"], [])
        self.assertEqual(inv["artifacts_with_forbidden_closure_true"], [])

    def test_decision_preserves_no_new_live_frontier(self):
        self.assertEqual(self.payload["status"], "P2707_POST_P2706_NO_NEW_LIVE_FRONTIER_CERTIFICATE")
        decision = self.payload["decision"]
        self.assertTrue(decision["no_new_live_frontier_certificate"])
        self.assertTrue(all(value is False for value in decision["negative_export_flags"].values()))
        self.assertIn("genuinely new strict typed object", decision["next_honest_step"])

    def test_outputs_and_docs_written(self):
        self.assertIn("P2707/S1657", MD.read_text(encoding="utf-8"))
        self.assertIn("P2707/S1657", (ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md").read_text(encoding="utf-8"))
        self.assertIn("P2707/S1657", (ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md").read_text(encoding="utf-8"))
        self.assertIn("P2707", (REPO / "AGENTS.md").read_text(encoding="utf-8"))


if __name__ == "__main__":
    unittest.main()
