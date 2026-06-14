from __future__ import annotations

import json
import subprocess
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
SCRIPT = ROOT / "p2710_s1660_finite_aut_z12_anti_inversion_orientation_character_source_test.py"
OUT = ROOT / "generated" / "p2710_s1660_finite_aut_z12_anti_inversion_orientation_character_source_test.json"
MD = ROOT / "generated" / "p2710_s1660_finite_aut_z12_anti_inversion_orientation_character_source_test.md"


class P2710AutZ12AntiInversionCharacterSourceTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        subprocess.run([sys.executable, str(SCRIPT)], cwd=ROOT, check=True)
        cls.payload = json.loads(OUT.read_text(encoding="utf-8"))

    def test_character_table_is_complete(self):
        summary = self.payload["character_summary"]
        self.assertEqual(summary["total_characters"], 4)
        self.assertEqual(summary["anti_inversion_characters"], 2)
        self.assertEqual(len(summary["anti_inversion_character_names"]), 2)
        self.assertTrue(all(row["homomorphism_valid"] for row in self.payload["character_table"]))

    def test_anti_inversion_characters_do_not_select_plus_omega(self):
        anti_rows = [row for row in self.payload["character_table"] if row["anti_inversion"]]
        self.assertEqual(len(anti_rows), 2)
        self.assertTrue(all(row["distinguishes_inversion_from_identity"] for row in anti_rows))
        self.assertTrue(all(not row["selects_plus_omega_without_source"] for row in anti_rows))

    def test_no_strict_source_exported(self):
        self.assertEqual(self.payload["status"], "P2710_AUT_Z12_ANTI_INVERSION_CHARACTER_EXISTS_BUT_NO_STRICT_SOURCE")
        self.assertTrue(all(not row["strict_anti_inversion_source_exported"] for row in self.payload["artifact_source_scan"]))
        self.assertTrue(all(not row["exports_orientation_source"] for row in self.payload["candidate_source_rows"]))
        decision = self.payload["decision"]
        self.assertTrue(decision["finite_anti_inversion_characters_exist"])
        self.assertFalse(decision["strict_orientation_character_source_exported"])
        self.assertTrue(all(value is False for value in decision["negative_export_flags"].values()))

    def test_outputs_and_docs_written(self):
        self.assertIn("P2710/S1660", MD.read_text(encoding="utf-8"))
        self.assertIn("P2710/S1660", (ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md").read_text(encoding="utf-8"))
        self.assertIn("P2710/S1660", (ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md").read_text(encoding="utf-8"))
        self.assertIn("P2710", (REPO / "AGENTS.md").read_text(encoding="utf-8"))


if __name__ == "__main__":
    unittest.main()
