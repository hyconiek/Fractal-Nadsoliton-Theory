from __future__ import annotations

import json
import subprocess
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent
SCRIPT = ROOT / "p2698_s1648_symmetry_breaking_direction_claim_reconciliation_audit.py"
OUT = ROOT / "generated" / "p2698_s1648_symmetry_breaking_direction_claim_reconciliation_audit.json"
MD = ROOT / "generated" / "p2698_s1648_symmetry_breaking_direction_claim_reconciliation_audit.md"


class P2698SymmetryBreakingDirectionClaimReconciliationTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        subprocess.run([sys.executable, str(SCRIPT)], cwd=ROOT, check=True)
        cls.payload = json.loads(OUT.read_text(encoding="utf-8"))

    def test_content_grep_covers_direction_claim_and_boundaries(self) -> None:
        audit = self.payload["content_grep"]
        self.assertIn("symmetry-breaking/direction", audit["mode"])
        for key in (
            "nadsoliton_information_ontology",
            "symmetry_breaking_direction_claims",
            "premise_based_scope_markers",
            "qw2191_open_markers",
            "p739_p740_nonexport",
            "no_false_promotion_boundaries",
        ):
            self.assertIn(key, audit["patterns"])
            self.assertGreater(audit["patterns"][key]["count"], 0)

    def test_state_reads_preserve_real_direction_but_keep_qw2191_open(self) -> None:
        state = self.payload["state_reads"]
        self.assertTrue(state["nadsoliton_ontology_respected"])
        self.assertTrue(state["directed_axis_orientation_present"])
        self.assertTrue(state["h37_sign_sensitive_directed_observable_exported"])
        self.assertTrue(state["h38_directed_continuation_selected"])
        self.assertTrue(state["h39_qw2191_still_open"])
        self.assertTrue(state["h41_qw2191_still_open"])
        self.assertTrue(state["p739_directed_lane_exported"])
        self.assertTrue(state["p739_lane_is_premise_based"])
        self.assertFalse(state["p739_pair12_upgrade_exported"])
        self.assertTrue(state["p740_sign_fixed_lane_exported"])
        self.assertTrue(state["p740_lane_is_gauge_only"])
        self.assertFalse(state["p740_pair12_upgrade_exported"])
        self.assertTrue(state["p1227_w1_witness_discharged_but_theory_open"])
        self.assertTrue(state["p1233_readiness_not_global_closure"])

    def test_reconciliation_matrix_blocks_silent_promotion(self) -> None:
        rows = self.payload["reconciliation_matrix"]
        self.assertEqual(len(rows), 4)
        self.assertTrue(all(row["support_is_real"] for row in rows))
        self.assertTrue(all(not row["promotion_allowed_now"] for row in rows))
        scopes = " ".join(row["scope"] for row in rows)
        self.assertIn("premise-based", scopes)
        self.assertIn("nonexport", scopes)

    def test_decision_preserves_no_new_live_frontier_certificate(self) -> None:
        decision = self.payload["decision"]
        self.assertTrue(decision["no_new_live_frontier_certificate_preserved"])
        self.assertIn("real direction/orientation", decision["acknowledgement"])
        self.assertIn("P2699", decision["next_honest_step"])
        self.assertTrue(all(value is False for value in self.payload["negative_export_flags"].values()))
        self.assertIn("P2698/S1648", MD.read_text(encoding="utf-8"))
        self.assertIn("P2698/S1648", (ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md").read_text(encoding="utf-8"))
        self.assertIn("P2698/S1648", (ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md").read_text(encoding="utf-8"))
        self.assertIn("P2698/S1648", (ROOT.parent / "AGENTS.md").read_text(encoding="utf-8"))


if __name__ == "__main__":
    unittest.main()
