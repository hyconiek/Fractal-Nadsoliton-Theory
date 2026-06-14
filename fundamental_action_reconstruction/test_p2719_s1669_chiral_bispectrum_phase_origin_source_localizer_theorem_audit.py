from __future__ import annotations

import json
import subprocess
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
SCRIPT = ROOT / "p2719_s1669_chiral_bispectrum_phase_origin_source_localizer_theorem_audit.py"
OUT = ROOT / "generated" / "p2719_s1669_chiral_bispectrum_phase_origin_source_localizer_theorem_audit.json"
MD = ROOT / "generated" / "p2719_s1669_chiral_bispectrum_phase_origin_source_localizer_theorem_audit.md"


class P2719ChiralBispectrumPhaseOriginSourceLocalizerTheoremAuditTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        subprocess.run([sys.executable, str(SCRIPT)], cwd=ROOT, check=True)
        cls.payload = json.loads(OUT.read_text(encoding="utf-8"))

    def test_phase_origin_theorem_remains_missing(self):
        self.assertEqual(self.payload["status"], "P2719_CHIRAL_BISPECTRUM_PHASE_ORIGIN_LOCALIZER_THEOREM_NO_UNLOCK")
        facts = self.payload["theorem_audit"]["facts"]
        self.assertTrue(facts["exact_formula_fixed"])
        self.assertTrue(facts["signed_marker_nonzero_and_orientation_separating"])
        self.assertFalse(facts["phase_origin_reference_internal_to_strict_artifacts"])
        self.assertFalse(facts["source_localizer_selects_one_origin_nonpremise"])
        self.assertFalse(facts["torsor_coupling_theorem_exported"])
        self.assertFalse(self.payload["theorem_audit"]["phase_origin_source_localizer_theorem_exported"])

    def test_no_closure_flags_are_exported(self):
        decision = self.payload["decision"]
        self.assertFalse(decision["phase_origin_source_localizer_theorem_exported"])
        self.assertFalse(decision["strict_chiral_bispectrum_source_exported"])
        self.assertFalse(decision["strict_mechanism_fixing_lambda_exported"])
        self.assertTrue(all(value is False for value in decision["negative_export_flags"].values()))

    def test_docs_written(self):
        self.assertIn("P2719/S1669", MD.read_text(encoding="utf-8"))
        self.assertIn("P2719/S1669", (ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md").read_text(encoding="utf-8"))
        self.assertIn("P2719/S1669", (ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md").read_text(encoding="utf-8"))
        self.assertIn("P2719", (REPO / "AGENTS.md").read_text(encoding="utf-8"))


if __name__ == "__main__":
    unittest.main()
