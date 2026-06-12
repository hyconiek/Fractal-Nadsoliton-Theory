from __future__ import annotations

import json
import subprocess
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent
SCRIPT = ROOT / "p2644_s1594_modified_compressed_inverse_hierarchy_successor_theorem.py"
OUT = ROOT / "generated" / "p2644_s1594_modified_compressed_inverse_hierarchy_successor_theorem.json"
MD = ROOT / "generated" / "p2644_s1594_modified_compressed_inverse_hierarchy_successor_theorem.md"


class P2644ModifiedCompressedInverseHierarchySuccessorTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        subprocess.run([sys.executable, str(SCRIPT)], cwd=ROOT, check=True)
        cls.payload = json.loads(OUT.read_text(encoding="utf-8"))

    def test_content_first_grep_is_present(self) -> None:
        audit = self.payload["semantic_rg_antiduplication_audit"]
        self.assertIn("content-first", audit["mode"])
        self.assertIn("modified_compressed_successor_content", audit["patterns"])
        self.assertIn("attention_suppression_content", audit["patterns"])
        self.assertTrue(all(data["count"] > 0 for data in audit["patterns"].values()))

    def test_suppression_factor_is_monotone_and_below_one(self) -> None:
        theorem = self.payload["compression_successor_theorem"]
        self.assertTrue(theorem["all_grid_suppression_below_one"])
        self.assertTrue(theorem["all_grid_derivative_negative"])
        self.assertAlmostEqual(theorem["d1_suppression"], 0.505, places=12)
        self.assertLess(theorem["d7_suppression"], theorem["d1_suppression"])
        self.assertLess(theorem["d7_over_d1_suppression_ratio"], 0.2)

    def test_role_successor_decision_is_modified_not_full_transfer(self) -> None:
        decision = self.payload["role_successor_decision"]
        self.assertTrue(decision["gates"]["distant_octave_suppressed_relative_to_near"])
        self.assertFalse(decision["gates"]["unchanged_inverse_hierarchy_reopened"])
        self.assertFalse(decision["gates"]["beta_source_available"])
        self.assertIn("UNCHANGED_INVERSE_HIERARCHY_REJECTED", decision["legacy_role_transfer_verdict"])
        self.assertFalse(decision["full_kernel_now"])
        self.assertFalse(decision["toe_closure_now"])
        self.assertTrue(all(value is False for value in self.payload["negative_export_flags"].values()))

    def test_docs_are_updated(self) -> None:
        self.assertIn("compression/locality-bias", self.payload["role_successor_decision"]["next_honest_step"])
        self.assertIn("P2644/S1594", MD.read_text(encoding="utf-8"))
        self.assertIn("P2644/S1594", (ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md").read_text(encoding="utf-8"))
        self.assertIn("P2644/S1594", (ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md").read_text(encoding="utf-8"))


if __name__ == "__main__":
    unittest.main()
