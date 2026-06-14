from __future__ import annotations

import json
import subprocess
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
SCRIPT = ROOT / "p2716_s1666_inversion_odd_pseudoscalar_source_acceptance_audit.py"
OUT = ROOT / "generated" / "p2716_s1666_inversion_odd_pseudoscalar_source_acceptance_audit.json"
MD = ROOT / "generated" / "p2716_s1666_inversion_odd_pseudoscalar_source_acceptance_audit.md"


class P2716InversionOddPseudoscalarSourceAcceptanceAuditTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        subprocess.run([sys.executable, str(SCRIPT)], cwd=ROOT, check=True)
        cls.payload = json.loads(OUT.read_text(encoding="utf-8"))

    def test_pseudoscalar_representation_has_two_equivariant_maps(self):
        self.assertEqual(self.payload["status"], "P2716_PSEUDOSCALAR_REPRESENTATION_ADMISSIBLE_BUT_NO_STRICT_SOURCE_VALUE")
        summary = self.payload["finite_summary"]
        self.assertEqual(summary["candidate_maps"], 4)
        self.assertEqual(summary["equivariant_maps"], 2)
        self.assertFalse(summary["strict_nonzero_pseudoscalar_source_exported"])
        self.assertEqual(len(self.payload["acceptance_rows"]), 2)
        self.assertTrue(all(row["representation_theoretically_admissible"] for row in self.payload["acceptance_rows"]))

    def test_no_source_value_means_no_lambda_fix_or_closure(self):
        self.assertTrue(all(not row["strict_nonzero_source_sign_exported"] for row in self.payload["acceptance_rows"]))
        self.assertTrue(all(not row["fixes_orientation_torsor_now"] for row in self.payload["acceptance_rows"]))
        decision = self.payload["decision"]
        self.assertTrue(decision["pseudoscalar_representation_can_couple_to_orientation_torsor"])
        self.assertFalse(decision["strict_pseudoscalar_or_chiral_source_exported"])
        self.assertFalse(decision["pseudoscalar_source_fixes_orientation_torsor"])
        self.assertTrue(all(value is False for value in decision["negative_export_flags"].values()))

    def test_docs_written(self):
        self.assertIn("P2716/S1666", MD.read_text(encoding="utf-8"))
        self.assertIn("P2716/S1666", (ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md").read_text(encoding="utf-8"))
        self.assertIn("P2716/S1666", (ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md").read_text(encoding="utf-8"))
        self.assertIn("P2716", (REPO / "AGENTS.md").read_text(encoding="utf-8"))


if __name__ == "__main__":
    unittest.main()
