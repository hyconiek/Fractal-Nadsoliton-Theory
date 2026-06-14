from __future__ import annotations

import json
import subprocess
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
SCRIPT = ROOT / "p2720_s1670_chiral_bispectrum_translation_orbit_phase_origin_localizer_no_go.py"
OUT = ROOT / "generated" / "p2720_s1670_chiral_bispectrum_translation_orbit_phase_origin_localizer_no_go.json"
MD = ROOT / "generated" / "p2720_s1670_chiral_bispectrum_translation_orbit_phase_origin_localizer_no_go.md"


class P2720ChiralBispectrumTranslationOrbitPhaseOriginLocalizerNoGoTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        subprocess.run([sys.executable, str(SCRIPT)], cwd=ROOT, check=True)
        cls.payload = json.loads(OUT.read_text(encoding="utf-8"))

    def test_translation_orbit_blocks_phase_origin_localizer(self):
        self.assertEqual(self.payload["status"], "P2720_TRANSLATION_ORBIT_PHASE_ORIGIN_LOCALIZER_NO_GO")
        audit = self.payload["translation_orbit_audit"]
        self.assertEqual(audit["accepted_localizer_count"], 0)
        self.assertTrue(audit["facts"]["exact_p2718_marker_rows_reused"])
        self.assertTrue(audit["facts"]["each_orientation_has_full_z12_source_orbit"])
        self.assertTrue(audit["facts"]["marker_constant_on_each_translation_orbit"])
        self.assertTrue(audit["facts"]["orientation_sign_still_recoverable"])
        self.assertFalse(audit["facts"]["source_origin_localizer_exported"])

    def test_orbit_rows_have_no_source_selector(self):
        for row in self.payload["translation_orbit_audit"]["orbit_rows"]:
            self.assertEqual(row["source_orbit_size"], 12)
            self.assertEqual(len(row["marker_values"]), 1)
            self.assertFalse(row["can_select_one_origin_without_external_label"])

    def test_no_closure_flags_are_exported(self):
        decision = self.payload["decision"]
        self.assertFalse(decision["translation_orbit_phase_origin_localizer_exported"])
        self.assertFalse(decision["nonpremise_source_origin_selected"])
        self.assertFalse(decision["strict_chiral_bispectrum_source_exported"])
        self.assertFalse(decision["strict_mechanism_fixing_lambda_exported"])
        self.assertTrue(all(value is False for value in decision["negative_export_flags"].values()))

    def test_docs_written(self):
        self.assertIn("P2720/S1670", MD.read_text(encoding="utf-8"))
        self.assertIn("P2720/S1670", (ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md").read_text(encoding="utf-8"))
        self.assertIn("P2720/S1670", (ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md").read_text(encoding="utf-8"))
        self.assertIn("P2720", (REPO / "AGENTS.md").read_text(encoding="utf-8"))


if __name__ == "__main__":
    unittest.main()
