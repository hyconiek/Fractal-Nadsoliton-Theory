from __future__ import annotations

import json
import subprocess
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
SCRIPT = ROOT / "p2718_s1668_chiral_bispectrum_explicit_signed_formula_torsor_coupling_audit.py"
OUT = ROOT / "generated" / "p2718_s1668_chiral_bispectrum_explicit_signed_formula_torsor_coupling_audit.json"
MD = ROOT / "generated" / "p2718_s1668_chiral_bispectrum_explicit_signed_formula_torsor_coupling_audit.md"


class P2718ChiralBispectrumExplicitSignedFormulaTorsorCouplingAuditTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        subprocess.run([sys.executable, str(SCRIPT)], cwd=ROOT, check=True)
        cls.payload = json.loads(OUT.read_text(encoding="utf-8"))

    def test_chiral_bispectrum_marker_is_nonzero_and_orientation_separating(self):
        self.assertEqual(self.payload["status"], "P2718_CHIRAL_BISPECTRUM_SIGNED_FORMULA_POSITIVE_BUT_NO_STRICT_TORSOR_SOURCE")
        summary = self.payload["finite_summary"]
        self.assertEqual(summary["checked_rows"], 24)
        self.assertTrue(summary["orientation_separating"])
        self.assertEqual(summary["orientation_marker_values"], {"-1": [2.0], "1": [-2.0]})
        self.assertTrue(all(row["nonzero_signed_marker"] for row in self.payload["marker_rows"]))

    def test_formula_is_not_accepted_without_phase_origin_and_torsor_coupling(self):
        acceptance = self.payload["acceptance_matrix"]
        self.assertTrue(acceptance["facts"]["explicit_signed_formula_computed"])
        self.assertTrue(acceptance["facts"]["nonzero_signed_value_on_all_rows"])
        self.assertFalse(acceptance["facts"]["phase_origin_reference_nonpremise"])
        self.assertFalse(acceptance["facts"]["coupling_to_p2708_p2714_torsor_exported"])
        self.assertFalse(acceptance["accepted_as_strict_pseudoscalar_source"])
        decision = self.payload["decision"]
        self.assertFalse(decision["chiral_bispectrum_accepted_as_strict_source"])
        self.assertFalse(decision["strict_mechanism_fixing_lambda_exported"])
        self.assertTrue(all(value is False for value in decision["negative_export_flags"].values()))

    def test_docs_written(self):
        self.assertIn("P2718/S1668", MD.read_text(encoding="utf-8"))
        self.assertIn("P2718/S1668", (ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md").read_text(encoding="utf-8"))
        self.assertIn("P2718/S1668", (ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md").read_text(encoding="utf-8"))
        self.assertIn("P2718", (REPO / "AGENTS.md").read_text(encoding="utf-8"))


if __name__ == "__main__":
    unittest.main()
