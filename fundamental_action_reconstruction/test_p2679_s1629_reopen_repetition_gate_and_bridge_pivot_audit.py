from __future__ import annotations

import json
import subprocess
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent
SCRIPT = ROOT / "p2679_s1629_reopen_repetition_gate_and_bridge_pivot_audit.py"
OUT = ROOT / "generated" / "p2679_s1629_reopen_repetition_gate_and_bridge_pivot_audit.json"
MD = ROOT / "generated" / "p2679_s1629_reopen_repetition_gate_and_bridge_pivot_audit.md"


class P2679ReopenRepetitionGateAndBridgePivotAuditTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        subprocess.run([sys.executable, str(SCRIPT)], cwd=ROOT, check=True)
        cls.payload = json.loads(OUT.read_text(encoding="utf-8"))

    def test_content_first_grep_covers_repeated_and_pivot_lanes(self) -> None:
        audit = self.payload["semantic_rg_antiduplication_audit"]
        self.assertIn("content-first", audit["mode"])
        for key in (
            "tau_src_pair12_repetition_content",
            "selector_orientation_repetition_content",
            "beta_tors_chi11_repetition_content",
            "legacy_strict_bridge_content",
            "damping_compression_content",
            "role_transfer_guard_content",
        ):
            self.assertIn(key, audit["patterns"])
            self.assertGreater(audit["patterns"][key]["count"], 0)

    def test_prior_research_ledger_blocks_repeats(self) -> None:
        ledger = {row["lane"]: row for row in self.payload["prior_research_ledger"]}
        self.assertIn("tau_src -> pair12 -> boundary-square", ledger)
        self.assertIn("selector/orientation source", ledger)
        self.assertIn("beta_tors -> chi11", ledger)
        self.assertIn("strict beta/damping source", ledger)
        self.assertTrue(all(not row["reopen_now"] for row in ledger.values()))
        self.assertIn("P2677", ledger["tau_src -> pair12 -> boundary-square"]["already_run"])
        self.assertIn("P2678", ledger["selector/orientation source"]["already_run"])
        self.assertIn("P2619", ledger["beta_tors -> chi11"]["already_run"])

    def test_reopen_gate_lattice_requires_new_evidence(self) -> None:
        lattice = self.payload["reopen_gate_lattice"]
        self.assertEqual(lattice["total_states"], 64)
        self.assertEqual(lattice["passing_states"], 1)
        self.assertEqual(lattice["hamming_distance_to_pass"], 2)
        self.assertFalse(lattice["may_reopen_old_lanes_now"])
        self.assertIn("new_object_after_p2678", lattice["missing_current_obligations"])
        self.assertIn("typed_precollapse_binding_or_bridge_source_exported", lattice["missing_current_obligations"])

    def test_bridge_pivot_matrix_selects_nonrepeat_lane(self) -> None:
        matrix = self.payload["bridge_pivot_matrix"]
        admissible = [row["candidate_next_lane"] for row in matrix if row["admissible_now"]]
        self.assertEqual(admissible, ["legacy->strict bridge source audit excluding selector/tau_src/beta_tors->chi11 repeats"])
        blocked = [row for row in matrix if not row["admissible_now"]]
        self.assertEqual(len(blocked), 3)
        self.assertTrue(any("tau_src" in row["candidate_next_lane"] for row in blocked))
        self.assertTrue(any("selector" in row["candidate_next_lane"] for row in blocked))
        self.assertTrue(any("beta_tors" in row["candidate_next_lane"] for row in blocked))

    def test_no_closure_and_docs_updated(self) -> None:
        decision = self.payload["closure_decision"]
        self.assertFalse(decision["old_lanes_reopened_now"])
        self.assertFalse(decision["new_reopening_evidence_after_p2678"])
        self.assertFalse(decision["s6_o3_reopened_now"])
        self.assertFalse(decision["o4_o5_allowed_now"])
        self.assertFalse(decision["role_bearing_ltotal_now"])
        self.assertFalse(decision["toe_closure_now"])
        self.assertTrue(all(value is False for value in self.payload["negative_export_flags"].values()))
        self.assertIn("P2679/S1629", MD.read_text(encoding="utf-8"))
        self.assertIn("P2679/S1629", (ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md").read_text(encoding="utf-8"))
        self.assertIn("P2679/S1629", (ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md").read_text(encoding="utf-8"))


if __name__ == "__main__":
    unittest.main()
