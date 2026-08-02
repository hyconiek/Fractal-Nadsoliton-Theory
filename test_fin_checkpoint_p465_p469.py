#!/usr/bin/env python3
"""Integrity tests for the FIN Release-10.43 checkpoint."""

import json
from pathlib import Path
import subprocess
import unittest


ROOT = Path(__file__).resolve().parent


class CheckpointP465P469Tests(unittest.TestCase):
    def test_required_artifacts(self) -> None:
        for name in (
            "FIN_Local_Research_Checkpoint_P465_P469_EN.md",
            "FIN_Local_Research_Checkpoint_P465_P469_EN.tex",
            "FIN_Local_Research_Checkpoint_P465_P469_EN.pdf",
            "FIN_Local_Research_Checkpoint_P465_P469_State.json",
            "FIN_Programs_465_468_469_Results.json",
            "FIN_Programs_465_468_469_P468_Rational_Dual.npz",
            "RELEASE_10_43_PROGRAMS_465_468_469.md",
            "FIN_Checkpoint_P465_P469_AGENTS_Guardrail.txt",
        ):
            self.assertTrue((ROOT / name).is_file(), name)

    def test_state_and_boundaries(self) -> None:
        state = json.loads(
            (ROOT / "FIN_Local_Research_Checkpoint_P465_P469_State.json").read_text(
                encoding="utf-8"
            )
        )
        self.assertEqual(state["next_selected_batch"], ["P471", "P472", "P473"])
        self.assertTrue(state["p465_full_simplex_unique"])
        self.assertLess(state["p468_global_gap"], 1e-8)
        self.assertFalse(state["p468_exact_optimizer_proved"])
        self.assertFalse(state["p469_exact_root_proved"])
        for key in (
            "physical_evidence_admitted",
            "selector_discharged",
            "dimensional_source_exported",
            "legacy_strict_bridge_complete",
            "legacy_role_transfer_started",
        ):
            self.assertFalse(state[key])

    def test_pdf(self) -> None:
        pdf = ROOT / "FIN_Local_Research_Checkpoint_P465_P469_EN.pdf"
        self.assertGreater(pdf.stat().st_size, 450_000)
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
            "0.52332810067104088",
            "41055151",
            "O176",
            "O177",
            "O178",
            "P471",
        ):
            self.assertIn(token, extracted)
        log = (ROOT / "FIN_Local_Research_Checkpoint_P465_P469_EN.log").read_text(
            encoding="utf-8", errors="replace"
        )
        self.assertNotIn("LaTeX Error", log)
        self.assertNotIn("Missing character", log)
        self.assertNotIn("Fatal error", log)

    def test_current_copy(self) -> None:
        self.assertEqual(
            (ROOT / "FIN_Local_Research_Checkpoint_P465_P469_EN.pdf").read_bytes(),
            (ROOT / "FIN_Current_Local_Research_Report_EN.pdf").read_bytes(),
        )


if __name__ == "__main__":
    unittest.main()
