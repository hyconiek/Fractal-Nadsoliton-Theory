#!/usr/bin/env python3
"""Integrity tests for the FIN Release-10.45 checkpoint."""

import json
from pathlib import Path
import subprocess
import unittest


ROOT = Path(__file__).resolve().parent


class CheckpointP474P479P483Tests(unittest.TestCase):
    def test_required_artifacts(self) -> None:
        for name in (
            "FIN_Local_Research_Checkpoint_P474_P479_P483_EN.md",
            "FIN_Local_Research_Checkpoint_P474_P479_P483_EN.tex",
            "FIN_Local_Research_Checkpoint_P474_P479_P483_EN.pdf",
            "FIN_Local_Research_Checkpoint_P474_P479_P483_State.json",
            "FIN_Programs_474_479_483_Results.json",
            "FIN_Programs_474_479_483_P474_Flat_Direction.npz",
            "FIN_Programs_474_479_483_P479_Riccati_Trace_Core.lean",
            "RELEASE_10_45_PROGRAMS_474_479_483.md",
            "FIN_Checkpoint_P474_P479_P483_AGENTS_Guardrail.txt",
        ):
            self.assertTrue((ROOT / name).is_file(), name)

    def test_state_and_boundaries(self) -> None:
        state = json.loads(
            (ROOT / "FIN_Local_Research_Checkpoint_P474_P479_P483_State.json").read_text(
                encoding="utf-8"
            )
        )
        self.assertEqual(state["next_batch"], ["P484", "P480", "P475"])
        self.assertEqual(state["programs"]["P474"]["status"],
                         "strong_evidence_exact_certificate_open")
        self.assertEqual(state["programs"]["P479"]["status"],
                         "formal_source_exported_local_toolchain_missing")
        self.assertEqual(state["programs"]["P483"]["status"],
                         "computer_assisted_proof")
        self.assertEqual(state["selector_QW_2191"], "open")
        self.assertEqual(state["dimensional_source"], "missing")
        self.assertFalse(state["network_used"])
        self.assertFalse(state["laboratory_data_used"])
        self.assertFalse(state["external_audit_used"])

    def test_pdf(self) -> None:
        pdf = ROOT / "FIN_Local_Research_Checkpoint_P474_P479_P483_EN.pdf"
        self.assertGreater(pdf.stat().st_size, 400_000)
        info = subprocess.run(
            ["pdfinfo", pdf.name], cwd=ROOT, capture_output=True, text=True, check=True
        ).stdout
        self.assertIn("Krzysztof Żuchowski", info)
        self.assertIn("A4", info)
        extracted = subprocess.run(
            ["pdftotext", "-layout", pdf.name, "-"],
            cwd=ROOT,
            capture_output=True,
            text=True,
            check=True,
        ).stdout
        for token in (
            "P474",
            "P479",
            "P483",
            "O182",
            "O183",
            "O184",
            "3 × 10−9",
            "P484",
            "QW-2191",
        ):
            self.assertIn(token, extracted)
        log = (ROOT / "FIN_Local_Research_Checkpoint_P474_P479_P483_EN.log").read_text(
            encoding="utf-8", errors="replace"
        )
        self.assertNotIn("LaTeX Error", log)
        self.assertNotIn("Missing character", log)
        self.assertNotIn("Fatal error", log)
        self.assertNotIn("Overfull", log)

    def test_current_copy(self) -> None:
        self.assertEqual(
            (ROOT / "FIN_Local_Research_Checkpoint_P474_P479_P483_EN.pdf").read_bytes(),
            (ROOT / "FIN_Current_Local_Research_Report_EN.pdf").read_bytes(),
        )


if __name__ == "__main__":
    unittest.main()
