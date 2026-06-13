from __future__ import annotations

import json
import subprocess
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent
SCRIPT = ROOT / "p2684_s1634_pair12_cycle_cut_semantic_invariant_provider_audit.py"
OUT = ROOT / "generated" / "p2684_s1634_pair12_cycle_cut_semantic_invariant_provider_audit.json"
MD = ROOT / "generated" / "p2684_s1634_pair12_cycle_cut_semantic_invariant_provider_audit.md"


class P2684Pair12CycleCutSemanticInvariantProviderAuditTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        subprocess.run([sys.executable, str(SCRIPT)], cwd=ROOT, check=True)
        cls.payload = json.loads(OUT.read_text(encoding="utf-8"))

    def test_content_grep_covers_cycle_cut_terms(self) -> None:
        audit = self.payload["content_grep"]
        self.assertIn("cycle-cut", audit["mode"])
        for key in (
            "cycle_cut_entry_point",
            "pre_collapse_guards",
            "semantic_invariant_provider_candidates",
            "known_near_miss_or_nonexport",
            "forbidden_closure_claims",
        ):
            self.assertIn(key, audit["patterns"])
            self.assertGreater(audit["patterns"][key]["count"], 0)

    def test_current_evidence_reads_p2683_and_p946_near_miss(self) -> None:
        evidence = self.payload["current_evidence"]
        self.assertTrue(evidence["p2683_loop_stop_present"])
        self.assertTrue(evidence["p946_no_current_lawful_bridge_supplier_found"])
        self.assertTrue(evidence["p946_near_miss_only"])
        self.assertTrue(evidence["f301_carrier_exists"])

    def test_obligation_matrix_blocks_cycle_cut_without_provider(self) -> None:
        matrix = self.payload["obligation_matrix"]
        self.assertFalse(matrix["semantic_invariant_provider_exported_now"])
        self.assertIn("chart_label_retaining_pair12_typed_seed_or_subinterface_actual_export", matrix["failed_obligations"])
        self.assertIn("nonconventional_semantic_invariant_or_provider", matrix["failed_obligations"])

    def test_policy_pivots_to_lagrangian_eom_reverse_closure(self) -> None:
        decision = self.payload["decision"]
        self.assertFalse(decision["semantic_invariant_provider_exported_now"])
        self.assertFalse(decision["lower_boundary_recursion_primary"])
        self.assertTrue(decision["strict_lagrangian_eom_reverse_closure_is_next"])
        self.assertFalse(decision["qw2191_discharged_now"])
        self.assertFalse(decision["toe_closed_now"])

    def test_docs_updated_without_closure_exports(self) -> None:
        self.assertTrue(all(value is False for value in self.payload["negative_export_flags"].values()))
        self.assertIn("P2684/S1634", MD.read_text(encoding="utf-8"))
        self.assertIn("P2684/S1634", (ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md").read_text(encoding="utf-8"))
        self.assertIn("P2684/S1634", (ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md").read_text(encoding="utf-8"))
        self.assertIn("P2684/S1634", (ROOT.parent / "AGENTS.md").read_text(encoding="utf-8"))


if __name__ == "__main__":
    unittest.main()
