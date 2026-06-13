from __future__ import annotations

import json
import subprocess
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent
SCRIPT = ROOT / "p2701_s1651_strict_sourced_symmetry_breaking_provider_inventory_matrix.py"
OUT = ROOT / "generated" / "p2701_s1651_strict_sourced_symmetry_breaking_provider_inventory_matrix.json"
MD = ROOT / "generated" / "p2701_s1651_strict_sourced_symmetry_breaking_provider_inventory_matrix.md"


class P2701StrictSourcedSymmetryBreakingProviderInventoryTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        subprocess.run([sys.executable, str(SCRIPT)], cwd=ROOT, check=True)
        cls.payload = json.loads(OUT.read_text(encoding="utf-8"))

    def test_content_grep_covers_provider_inventory_boundary(self) -> None:
        audit = self.payload["content_grep"]
        self.assertIn("P2701", audit["mode"])
        for key in ("p2700_boundary", "provider_search_terms", "known_premise_boundaries", "forbidden_promotions"):
            self.assertIn(key, audit["patterns"])
            self.assertGreater(audit["patterns"][key]["count"], 0)

    def test_state_reads_preserve_prior_no_go_boundaries(self) -> None:
        state = self.payload["state_reads"]
        self.assertTrue(state["p2700_bounded_no_go"])
        self.assertTrue(state["p2700_forbids_replay"])
        self.assertTrue(state["p2697_no_new_live_frontier_certificate"])
        self.assertTrue(state["h39_qw2191_still_open"])
        self.assertTrue(state["p739_pair12_upgrade_unexported"])
        self.assertTrue(state["p740_pair12_upgrade_unexported"])

    def test_inventory_scans_candidates_but_accepts_no_provider(self) -> None:
        inv = self.payload["provider_inventory"]
        self.assertGreater(inv["generated_json_files_scanned"], 100)
        self.assertGreater(inv["candidate_files_with_selector_direction_terms"], 0)
        self.assertEqual(inv["accepted_provider_count"], 0)
        self.assertEqual(inv["accepted_providers"], [])
        self.assertGreater(inv["blocked_candidate_count"], 0)
        self.assertGreater(len(inv["top_blocked_candidates"]), 0)

    def test_acceptance_matrix_and_decision_preserve_no_false_pass(self) -> None:
        rows = self.payload["acceptance_matrix"]
        self.assertEqual(len(rows), 4)
        self.assertTrue(all(not row["passes"] for row in rows))
        decision = self.payload["decision"]
        self.assertTrue(decision["bounded_no_go_now"])
        self.assertIn("no-new-live-frontier", decision["next_honest_step"])
        self.assertTrue(all(value is False for value in self.payload["negative_export_flags"].values()))
        self.assertIn("P2701/S1651", MD.read_text(encoding="utf-8"))
        self.assertIn("P2701/S1651", (ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md").read_text(encoding="utf-8"))
        self.assertIn("P2701/S1651", (ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md").read_text(encoding="utf-8"))
        self.assertIn("P2701/S1651", (ROOT.parent / "AGENTS.md").read_text(encoding="utf-8"))


if __name__ == "__main__":
    unittest.main()
