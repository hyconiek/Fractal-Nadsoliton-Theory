#!/usr/bin/env python3
"""Integrity tests for the FIN Release-10.38 durable checkpoint."""

from __future__ import annotations

import json
from pathlib import Path
import subprocess
import unittest


ROOT = Path(__file__).resolve().parent


class CheckpointP445P447Tests(unittest.TestCase):
    def test_required_artifacts_exist(self) -> None:
        for name in (
            "FIN_Local_Research_Checkpoint_P445_P447_EN.md",
            "FIN_Local_Research_Checkpoint_P445_P447_EN.tex",
            "FIN_Local_Research_Checkpoint_P445_P447_EN.pdf",
            "FIN_Current_Local_Research_Report_EN.pdf",
            "FIN_Local_Research_Checkpoint_P445_P447_State.json",
            "FIN_Programs_445_447_Results.json",
            "RELEASE_10_38_PROGRAMS_445_447.md",
            "FIN_Checkpoint_P445_P447_AGENTS_Guardrail.txt",
        ):
            self.assertTrue((ROOT / name).is_file(), name)

    def test_state_preserves_boundaries(self) -> None:
        state = json.loads(
            (ROOT / "FIN_Local_Research_Checkpoint_P445_P447_State.json").read_text(
                encoding="utf-8"
            )
        )
        self.assertEqual(state["latest_completed_program"], "P447")
        self.assertEqual(state["next_selected_batch"], ["P448", "P449", "P450"])
        self.assertTrue(state["p445_general_two_slot_tester_reduction_proved"])
        self.assertTrue(state["p445_arbitrary_intermediate_memory_included"])
        self.assertFalse(state["p445_three_or_more_slots_proved"])
        self.assertTrue(state["p446_palindromic_global_certificate_proved"])
        self.assertLess(state["p446_palindromic_gap"], 1e-3)
        self.assertFalse(state["p446_full_simplex_globality_proved"])
        self.assertTrue(state["p447_complete_p429_boxes_propagated"])
        self.assertTrue(state["p447_strict_mse_order_certified"])
        self.assertFalse(state["p447_detector_calibration_supplied"])
        self.assertFalse(state["physical_evidence_admitted"])
        self.assertFalse(state["selector_discharged"])
        self.assertFalse(state["dimensional_source_exported"])
        self.assertFalse(state["legacy_strict_bridge_complete"])

    def test_report_contains_required_boundaries(self) -> None:
        text = (ROOT / "FIN_Local_Research_Checkpoint_P445_P447_EN.md").read_text(
            encoding="utf-8"
        )
        for token in (
            "O161 -- Diagonal-Support Comb Reduction",
            "O162 -- Palindromic Heralded-Code Global Upper Certificate",
            "O163 -- Certified Detector-Allocation Tube",
            "8\\sqrt2}{25}",
            "9.924652413310642\\times10^{-4}",
            "QW-2191 is therefore unchanged",
            "Full-simplex globality remains open",
            "P448\\to P449\\to P450",
            "Restart instructions",
        ):
            self.assertIn(token, text)

    def test_pdf_metadata_text_figure_and_clean_log(self) -> None:
        pdf = ROOT / "FIN_Local_Research_Checkpoint_P445_P447_EN.pdf"
        self.assertGreater(pdf.stat().st_size, 450_000)
        info = subprocess.run(
            ["pdfinfo", pdf.name], cwd=ROOT, capture_output=True, text=True,
            timeout=30, check=True,
        ).stdout
        self.assertIn("FIN Local Research Checkpoint P445", info)
        self.assertIn("Krzysztof Żuchowski", info)
        pages = next(
            int(line.split(":", 1)[1].strip())
            for line in info.splitlines() if line.startswith("Pages:")
        )
        self.assertGreaterEqual(pages, 25)

        extracted = subprocess.run(
            ["pdftotext", "-layout", pdf.name, "-"], cwd=ROOT,
            capture_output=True, text=True, timeout=30, check=True,
        ).stdout
        for token in (
            "Krzysztof Żuchowski",
            "0.4525483399593904",
            "0.4642707484333634",
            "3.2548\\times10",
            "P 448 → P 449 → P 450",
        ):
            self.assertIn(token, extracted)
        self.assertGreater(
            (ROOT / "FIN_Programs_445_447_Figures/p445_p447_exact_reduction_and_intervals.png").stat().st_size,
            100_000,
        )
        log = (ROOT / "FIN_Local_Research_Checkpoint_P445_P447_EN.log").read_text(
            encoding="utf-8", errors="replace"
        )
        self.assertNotIn("Missing character", log)
        self.assertNotIn("LaTeX Error", log)

    def test_current_report_matches_latest_checkpoint(self) -> None:
        candidates = (
            ROOT / "FIN_Local_Research_Checkpoint_P458_P464_EN.pdf",
            ROOT / "FIN_Local_Research_Checkpoint_P454_P457_EN.pdf",
            ROOT / "FIN_Local_Research_Checkpoint_P451_P453_EN.pdf",
            ROOT / "FIN_Local_Research_Checkpoint_P448_P450_EN.pdf",
            ROOT / "FIN_Local_Research_Checkpoint_P445_P447_EN.pdf",
        )
        checkpoint_path = next(path for path in candidates if path.is_file())
        checkpoint = checkpoint_path.read_bytes()
        current = (ROOT / "FIN_Current_Local_Research_Report_EN.pdf").read_bytes()
        self.assertEqual(checkpoint, current)


if __name__ == "__main__":
    unittest.main()
