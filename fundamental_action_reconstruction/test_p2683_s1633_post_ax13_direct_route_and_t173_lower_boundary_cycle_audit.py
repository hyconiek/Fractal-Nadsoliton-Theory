from __future__ import annotations

import json
import subprocess
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent
SCRIPT = ROOT / "p2683_s1633_post_ax13_direct_route_and_t173_lower_boundary_cycle_audit.py"
OUT = ROOT / "generated" / "p2683_s1633_post_ax13_direct_route_and_t173_lower_boundary_cycle_audit.json"
MD = ROOT / "generated" / "p2683_s1633_post_ax13_direct_route_and_t173_lower_boundary_cycle_audit.md"


class P2683PostAx13DirectRouteAndT173LowerBoundaryCycleAuditTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        subprocess.run([sys.executable, str(SCRIPT)], cwd=ROOT, check=True)
        cls.payload = json.loads(OUT.read_text(encoding="utf-8"))

    def test_content_grep_covers_current_frontiers(self) -> None:
        audit = self.payload["content_grep"]
        self.assertIn("content-first", audit["mode"])
        for key in (
            "ax13_p51_direct_route_after_target_eom",
            "direct_route_negative_freeze",
            "h37_t171_premise_based_directed_state",
            "strict_core_upgrade_nonexport",
            "lower_boundary_recursion",
            "forbidden_closure_claims",
        ):
            self.assertIn(key, audit["patterns"])
            self.assertGreater(audit["patterns"][key]["count"], 0)

    def test_current_state_corrects_p2682_target(self) -> None:
        state = self.payload["current_artifact_state"]
        self.assertTrue(state["p2682_exists"])
        self.assertTrue(state["ax13_target_eom_closed_external"])
        self.assertFalse(state["ax13_strict_core_promotion"])
        self.assertFalse(state["p51_full_closure_pass"])
        self.assertTrue(state["p631_direct_route_negative_freeze_selected"])

    def test_h37_exists_but_pair12_strict_core_upgrade_still_blocked(self) -> None:
        state = self.payload["current_artifact_state"]
        self.assertTrue(state["h37_premise_based_directed_state_exported"])
        self.assertFalse(state["p739_directed_state_strict_core_upgrade_exported"])
        self.assertFalse(state["p740_sign_fixed_closure_strict_core_upgrade_exported"])
        self.assertFalse(state["p739_pair12_split_upgrades"])
        self.assertFalse(state["p740_pair12_split_upgrades"])

    def test_lower_boundary_cycle_detected_and_not_primary_without_invariant(self) -> None:
        sequence = self.payload["lower_boundary_sequence"]
        self.assertGreaterEqual(len(sequence), 8)
        lattice = self.payload["recursion_lattice"]
        self.assertTrue(lattice["lower_boundary_loop_detected"])
        self.assertFalse(lattice["continue_lower_boundary_without_new_invariant_admissible"])
        self.assertFalse(lattice["current_state"]["new_semantic_invariant_exported"])
        self.assertFalse(lattice["current_state"]["continue_lower_boundary_as_primary"])

    def test_no_closure_and_docs_updated(self) -> None:
        decision = self.payload["closure_decision"]
        self.assertFalse(decision["target_eom_work_reopened"])
        self.assertFalse(decision["h37_t171_replayed_as_missing"])
        self.assertFalse(decision["lower_boundary_recursion_primary"])
        self.assertFalse(decision["toe_closed_now"])
        self.assertTrue(all(value is False for value in self.payload["negative_export_flags"].values()))
        self.assertIn("P2683/S1633", MD.read_text(encoding="utf-8"))
        self.assertIn("P2683/S1633", (ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md").read_text(encoding="utf-8"))
        self.assertIn("P2683/S1633", (ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md").read_text(encoding="utf-8"))


if __name__ == "__main__":
    unittest.main()
