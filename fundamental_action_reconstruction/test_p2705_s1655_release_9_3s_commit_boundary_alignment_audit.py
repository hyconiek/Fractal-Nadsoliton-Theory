from __future__ import annotations

import json
import subprocess
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
SCRIPT = ROOT / "p2705_s1655_release_9_3s_commit_boundary_alignment_audit.py"
OUT = ROOT / "generated" / "p2705_s1655_release_9_3s_commit_boundary_alignment_audit.json"
MD = ROOT / "generated" / "p2705_s1655_release_9_3s_commit_boundary_alignment_audit.md"
COMMIT = "8d48faf012f87721d01a692fd7e3888461d4e6d2"


class P2705Release93sCommitBoundaryAuditTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        subprocess.run([sys.executable, str(SCRIPT)], cwd=ROOT, check=True)
        cls.payload = json.loads(OUT.read_text(encoding="utf-8"))

    def test_commit_identity_and_scope(self):
        meta = self.payload["commit_metadata"]
        self.assertEqual(meta["commit"], COMMIT)
        self.assertIn("P2377/P2378", meta["message"])
        self.assertTrue(any("p2377" in path for path in meta["changed_files"]))
        self.assertTrue(any("p2378" in path for path in meta["changed_files"]))

    def test_alignment_matrix_passes_without_selector_unlock(self):
        self.assertTrue(all(row["passes"] for row in self.payload["alignment_matrix"]))
        decision = self.payload["decision"]
        self.assertFalse(decision["release_9_3s_pointer_unblocks_current_stage"])
        self.assertTrue(all(value is False for value in decision["negative_export_flags"].values()))
        self.assertIn("P2706", decision["next_honest_step"])

    def test_numeric_boundary_records_insufficiency(self):
        numeric = self.payload["p2377_p2378_numeric_boundary"]
        self.assertTrue(numeric["p2378_unit_mass_insufficient_on_rectangle"])
        self.assertIn("not derived", numeric["p2377_normalization_status"])

    def test_outputs_and_docs_written(self):
        self.assertEqual(self.payload["status"], "P2705_RELEASE_9_3S_POINTER_AUDITED_NO_CURRENT_UNLOCK")
        self.assertIn("P2705/S1655", MD.read_text(encoding="utf-8"))
        self.assertIn("P2705/S1655", (ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md").read_text(encoding="utf-8"))
        self.assertIn("P2705/S1655", (ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md").read_text(encoding="utf-8"))
        self.assertIn("P2705", (REPO / "AGENTS.md").read_text(encoding="utf-8"))


if __name__ == "__main__":
    unittest.main()
