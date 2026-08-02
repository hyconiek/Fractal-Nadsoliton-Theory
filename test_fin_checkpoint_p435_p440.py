#!/usr/bin/env python3
"""Integrity tests for the FIN Release-10.37 durable checkpoint."""

from __future__ import annotations

import json
from pathlib import Path
import subprocess
import unittest


ROOT = Path(__file__).resolve().parent


class CheckpointP435P440Tests(unittest.TestCase):
    def test_required_artifacts_exist(self) -> None:
        for name in (
            "FIN_Local_Research_Checkpoint_P435_P440_EN.md",
            "FIN_Local_Research_Checkpoint_P435_P440_EN.tex",
            "FIN_Local_Research_Checkpoint_P435_P440_EN.pdf",
            "FIN_Current_Local_Research_Report_EN.pdf",
            "FIN_Local_Research_Checkpoint_P435_P440_State.json",
            "FIN_Programs_435_436_440_Results.json",
            "RELEASE_10_37_PROGRAMS_435_436_440.md",
            "FIN_Checkpoint_P435_P440_AGENTS_Guardrail.txt",
        ):
            self.assertTrue((ROOT / name).is_file(), name)

    def test_state_preserves_boundaries(self) -> None:
        state = json.loads(
            (ROOT / "FIN_Local_Research_Checkpoint_P435_P440_State.json").read_text(
                encoding="utf-8"
            )
        )
        self.assertEqual(state["latest_completed_program"], "P440")
        self.assertEqual(state["next_selected_batch"], ["P445", "P446", "P447"])
        self.assertTrue(state["p435_one_slot_sdp_exact"])
        self.assertFalse(state["p435_two_slot_matching_dual_exported"])
        self.assertTrue(state["p436_positive_rational_code_gain_certified"])
        self.assertGreater(state["p436_certified_gain_lower"], 0.0225)
        self.assertTrue(state["p440_conditional_minimax_law_proved"])
        self.assertFalse(state["p440_detector_calibration_supplied"])
        self.assertFalse(state["physical_evidence_admitted"])
        self.assertFalse(state["selector_discharged"])
        self.assertFalse(state["dimensional_source_exported"])
        self.assertFalse(state["legacy_strict_bridge_complete"])

    def test_report_contains_required_boundaries(self) -> None:
        text = (ROOT / "FIN_Local_Research_Checkpoint_P435_P440_EN.md").read_text(
            encoding="utf-8"
        )
        for token in (
            "O158 -- Process-Tester SDP Admission Certificate",
            "O159 -- Rational Heralded-Code Gain Certificate",
            "O160 -- Detector-Envelope Minimax Jordan Sampler",
            "0.022572776021405654",
            "1,909,937",
            "QW-2191 remains open",
            "No complete legacy-to-strict bridge",
            "P445\\to P446\\to P447",
            "Restart instructions",
        ):
            self.assertIn(token, text)

    def test_pdf_metadata_text_and_figure(self) -> None:
        pdf = ROOT / "FIN_Local_Research_Checkpoint_P435_P440_EN.pdf"
        self.assertGreater(pdf.stat().st_size, 450_000)
        info = subprocess.run(
            ["pdfinfo", pdf.name], cwd=ROOT, capture_output=True, text=True,
            timeout=30, check=True,
        ).stdout
        self.assertIn("FIN Local Research Checkpoint P435", info)
        self.assertIn("Krzysztof Żuchowski", info)
        pages = next(
            int(line.split(":", 1)[1].strip())
            for line in info.splitlines() if line.startswith("Pages:")
        )
        self.assertGreaterEqual(pages, 24)

        extracted = subprocess.run(
            ["pdftotext", "-layout", pdf.name, "-"], cwd=ROOT,
            capture_output=True, text=True, timeout=30, check=True,
        ).stdout
        for token in (
            "Krzysztof Żuchowski",
            "0.022572776021405654",
            "1,909,937",
            "genuine two-slot comb",
            "P445/P446/P447",
        ):
            self.assertIn(token, extracted)
        self.assertGreater(
            (ROOT / "FIN_Programs_435_436_440_Figures/p435_p436_p440_certificates.png").stat().st_size,
            100_000,
        )

    def test_current_report_matches_checkpoint(self) -> None:
        candidates = (
            ROOT / "FIN_Local_Research_Checkpoint_P458_P464_EN.pdf",
            ROOT / "FIN_Local_Research_Checkpoint_P454_P457_EN.pdf",
            ROOT / "FIN_Local_Research_Checkpoint_P451_P453_EN.pdf",
            ROOT / "FIN_Local_Research_Checkpoint_P448_P450_EN.pdf",
            ROOT / "FIN_Local_Research_Checkpoint_P445_P447_EN.pdf",
            ROOT / "FIN_Local_Research_Checkpoint_P435_P440_EN.pdf",
        )
        checkpoint_path = next(path for path in candidates if path.is_file())
        checkpoint = checkpoint_path.read_bytes()
        current = (ROOT / "FIN_Current_Local_Research_Report_EN.pdf").read_bytes()
        self.assertEqual(checkpoint, current)


if __name__ == "__main__":
    unittest.main()
