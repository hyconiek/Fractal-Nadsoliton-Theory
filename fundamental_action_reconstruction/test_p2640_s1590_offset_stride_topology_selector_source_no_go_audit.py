from __future__ import annotations

import json
import subprocess
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent
SCRIPT = ROOT / "p2640_s1590_offset_stride_topology_selector_source_no_go_audit.py"
OUT = ROOT / "generated" / "p2640_s1590_offset_stride_topology_selector_source_no_go_audit.json"
MD = ROOT / "generated" / "p2640_s1590_offset_stride_topology_selector_source_no_go_audit.md"


class P2640OffsetStrideTopologySelectorSourceNoGoAuditTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        subprocess.run([sys.executable, str(SCRIPT)], cwd=ROOT, check=True)
        cls.payload = json.loads(OUT.read_text(encoding="utf-8"))

    def test_latest_commit_and_content_first_grep_are_present(self) -> None:
        self.assertGreaterEqual(len(self.payload["latest_commit_audit"]), 3)
        audit = self.payload["semantic_rg_antiduplication_audit"]
        self.assertIn("content-first", audit["mode"])
        self.assertIn("offset_stride_source_content", audit["patterns"])
        self.assertIn("topology_selector_carrier_content", audit["patterns"])
        self.assertTrue(all(data["count"] > 0 for data in audit["patterns"].values()))

    def test_topology_selector_sources_do_not_export_canonical_pair(self) -> None:
        ledger = self.payload["topology_selector_source_ledger"]
        self.assertGreaterEqual(len(ledger), 6)
        self.assertTrue(any(row["id"] == "phase12_aut_z12_quotient_parity" for row in ledger))
        self.assertTrue(any(row["id"] == "mode_index_uniqueness_obstruction" for row in ledger))
        self.assertFalse(any(row["strict_offset_stride_source_pass"] for row in ledger))

    def test_p2639_role_like_candidates_remain_unsourced(self) -> None:
        role_like = self.payload["p2639_role_like_candidates"]
        pairs = {(row["k0"], row["stride"]) for row in role_like}
        self.assertIn((4, 3), pairs)
        self.assertIn((10, 6), pairs)
        matrix = self.payload["compatibility_matrix"]
        self.assertFalse(any(row["licensed_without_hidden_selector"] for row in matrix))
        self.assertTrue(any(row["has_numeric_stride_hint"] for row in matrix))

    def test_closure_decision_does_not_promote(self) -> None:
        decision = self.payload["closure_decision"]
        self.assertTrue(decision["gates"]["p2639_role_like_lifts_exist"])
        self.assertTrue(decision["gates"]["repo_exports_topology_selector_carriers"])
        self.assertFalse(decision["gates"]["some_source_exports_canonical_k0_and_stride"])
        self.assertFalse(decision["promote_offset_stride_to_bridge_completion"])
        self.assertFalse(decision["full_kernel_now"])
        self.assertTrue(all(value is False for value in self.payload["negative_export_flags"].values()))

    def test_markdown_and_docs_updated(self) -> None:
        self.assertIn("P2640/S1590", MD.read_text(encoding="utf-8"))
        self.assertIn("P2640/S1590", (ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md").read_text(encoding="utf-8"))
        self.assertIn("P2640/S1590", (ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md").read_text(encoding="utf-8"))
        self.assertIn("quotient-safe successor", self.payload["closure_decision"]["next_honest_step"])


if __name__ == "__main__":
    unittest.main()
