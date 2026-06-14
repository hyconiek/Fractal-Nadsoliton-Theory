from __future__ import annotations

import json
import subprocess
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
SCRIPT = ROOT / "p2724_s1674_post_p2723_commit_intake_next_honest_step_certificate.py"
OUT = ROOT / "generated" / "p2724_s1674_post_p2723_commit_intake_next_honest_step_certificate.json"
MD = ROOT / "generated" / "p2724_s1674_post_p2723_commit_intake_next_honest_step_certificate.md"


class P2724PostP2723CommitIntakeNextHonestStepCertificateTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        subprocess.run([sys.executable, str(SCRIPT)], cwd=ROOT, check=True)
        cls.payload = json.loads(OUT.read_text(encoding="utf-8"))

    def test_commit_intake_preserves_no_new_frontier(self):
        self.assertEqual(cls_status := self.payload["status"], "P2724_POST_P2723_COMMIT_INTAKE_NO_NEW_LIVE_FRONTIER_CERTIFICATE")
        self.assertIn("88eb860b1658ac5b648253fa65dd83bd4abbe922", self.payload["commit"])
        self.assertTrue(self.payload["commit_summary"]["available"])
        self.assertTrue(self.payload["commit_summary"]["mentions_p2719_to_p2723"], cls_status)
        self.assertFalse(self.payload["intake"]["new_strict_dynamic_chiral_source_artifact_supplied_by_commit_intake"])

    def test_required_next_object_is_precise(self):
        required = self.payload["intake"]["required_new_object"]
        self.assertEqual(required["name"], "strict_dynamic_chiral_source_artifact_coupled_to_P2721")
        self.assertEqual(len(required["requirements"]), 5)
        self.assertTrue(any("computable nonzero signed value" in item for item in required["requirements"]))
        self.assertTrue(any("P2721" in item for item in required["requirements"]))

    def test_no_closure_flags_are_exported(self):
        flags = self.payload["decision"]["negative_export_flags"]
        self.assertTrue(flags)
        self.assertTrue(all(value is False for value in flags.values()))
        self.assertIn("L_total promotion", self.payload["intake"]["closed_or_blocked_lanes"])
        self.assertIn("ToE closure", self.payload["intake"]["closed_or_blocked_lanes"])

    def test_docs_written(self):
        self.assertIn("P2724/S1674", MD.read_text(encoding="utf-8"))
        self.assertIn("P2724/S1674", (ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md").read_text(encoding="utf-8"))
        self.assertIn("P2724/S1674", (ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md").read_text(encoding="utf-8"))
        self.assertIn("P2724", (REPO / "AGENTS.md").read_text(encoding="utf-8"))


if __name__ == "__main__":
    unittest.main()
