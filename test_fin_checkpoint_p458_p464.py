#!/usr/bin/env python3
"""Integrity tests for the FIN Release-10.42 checkpoint."""

import json
from pathlib import Path
import subprocess
import unittest

ROOT = Path(__file__).resolve().parent


class CheckpointP458P464Tests(unittest.TestCase):
    def test_required_artifacts(self) -> None:
        for name in (
            "FIN_Local_Research_Checkpoint_P458_P464_EN.md",
            "FIN_Local_Research_Checkpoint_P458_P464_EN.tex",
            "FIN_Local_Research_Checkpoint_P458_P464_EN.pdf",
            "FIN_Local_Research_Checkpoint_P458_P464_State.json",
            "FIN_Programs_458_459_464_Results.json",
            "FIN_Programs_458_459_464_P464_Rational_Dual.npz",
            "RELEASE_10_42_PROGRAMS_458_459_464.md",
            "FIN_Checkpoint_P458_P464_AGENTS_Guardrail.txt",
        ):
            self.assertTrue((ROOT/name).is_file(), name)

    def test_state_and_report_boundaries(self) -> None:
        state = json.loads((ROOT/"FIN_Local_Research_Checkpoint_P458_P464_State.json").read_text(encoding="utf-8"))
        self.assertEqual(state["next_selected_batch"], ["P465", "P468", "P469"])
        self.assertLess(state["p464_global_gap"], 2.3e-7)
        self.assertTrue(state["p458_palindromic_unique"])
        self.assertFalse(state["p458_full_simplex_unique"])
        for key in ("physical_evidence_admitted", "selector_discharged", "dimensional_source_exported", "legacy_strict_bridge_complete"):
            self.assertFalse(state[key])

    def test_pdf(self) -> None:
        pdf = ROOT/"FIN_Local_Research_Checkpoint_P458_P464_EN.pdf"
        self.assertGreater(pdf.stat().st_size, 370_000)
        info = subprocess.run(["pdfinfo",pdf.name], cwd=ROOT, capture_output=True, text=True, check=True).stdout
        self.assertIn("Krzysztof Żuchowski", info)
        self.assertIn("A4", info)
        extracted = subprocess.run(["pdftotext","-layout",pdf.name,"-"], cwd=ROOT, capture_output=True, text=True, check=True).stdout
        for token in ("0.52332832027103", "O173", "O174", "O175", "P465"):
            self.assertIn(token, extracted)
        log = (ROOT/"FIN_Local_Research_Checkpoint_P458_P464_EN.log").read_text(encoding="utf-8",errors="replace")
        self.assertNotIn("LaTeX Error", log)
        self.assertNotIn("Missing character", log)

    def test_current_copy(self) -> None:
        self.assertEqual((ROOT/"FIN_Local_Research_Checkpoint_P458_P464_EN.pdf").read_bytes(), (ROOT/"FIN_Current_Local_Research_Report_EN.pdf").read_bytes())


if __name__ == "__main__":
    unittest.main()
