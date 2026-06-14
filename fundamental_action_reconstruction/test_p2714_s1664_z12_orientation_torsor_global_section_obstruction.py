from __future__ import annotations

import json
import subprocess
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
SCRIPT = ROOT / "p2714_s1664_z12_orientation_torsor_global_section_obstruction.py"
OUT = ROOT / "generated" / "p2714_s1664_z12_orientation_torsor_global_section_obstruction.json"
MD = ROOT / "generated" / "p2714_s1664_z12_orientation_torsor_global_section_obstruction.md"


class P2714Z12OrientationTorsorGlobalSectionObstructionTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        subprocess.run([sys.executable, str(SCRIPT)], cwd=ROOT, check=True)
        cls.payload = json.loads(OUT.read_text(encoding="utf-8"))

    def test_aut_action_has_orientation_reversers(self):
        self.assertEqual(self.payload["status"], "P2714_ORIENTATION_TORSOR_GLOBAL_SECTION_OBSTRUCTION_NO_STRICT_LAMBDA_FIX")
        actions = {row["unit"]: row for row in self.payload["action_table"]}
        self.assertFalse(actions[1]["orientation_reversing"])
        self.assertFalse(actions[5]["orientation_reversing"])
        self.assertTrue(actions[7]["orientation_reversing"])
        self.assertTrue(actions[11]["orientation_reversing"])
        self.assertEqual(actions[11]["image_of_plus_omega"], "-omega")

    def test_no_aut_compatible_global_section_or_closure_flags(self):
        summary = self.payload["torsor_summary"]
        self.assertEqual(summary["global_section_count"], 0)
        self.assertFalse(summary["has_aut_compatible_global_section"])
        self.assertTrue(all(not row["aut_compatible"] for row in self.payload["global_section_candidates"]))
        decision = self.payload["decision"]
        self.assertTrue(decision["new_typed_candidate_tested"])
        self.assertFalse(decision["orientation_torsor_global_section_exported"])
        self.assertFalse(decision["strict_mechanism_fixing_lambda_exported"])
        self.assertTrue(all(value is False for value in decision["negative_export_flags"].values()))

    def test_prior_boundary_and_docs_written(self):
        self.assertTrue(self.payload["prior_boundary_check"]["p2713_no_new_frontier_preserved_before_this_candidate"])
        self.assertIn("P2714/S1664", MD.read_text(encoding="utf-8"))
        self.assertIn("P2714/S1664", (ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md").read_text(encoding="utf-8"))
        self.assertIn("P2714/S1664", (ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md").read_text(encoding="utf-8"))
        self.assertIn("P2714", (REPO / "AGENTS.md").read_text(encoding="utf-8"))


if __name__ == "__main__":
    unittest.main()
