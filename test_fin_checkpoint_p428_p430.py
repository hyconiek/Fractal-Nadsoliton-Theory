#!/usr/bin/env python3
"""Integrity tests for the FIN P428--P430 durable PDF checkpoint."""

from __future__ import annotations

import json
from pathlib import Path
import subprocess
import unittest


ROOT = Path(__file__).resolve().parent


class CheckpointP428P430Tests(unittest.TestCase):
    def test_required_checkpoint_artifacts_exist(self) -> None:
        for name in (
            "FIN_Local_Research_Checkpoint_P428_P430_EN.md",
            "FIN_Local_Research_Checkpoint_P428_P430_EN.tex",
            "FIN_Local_Research_Checkpoint_P428_P430_EN.pdf",
            "FIN_Current_Local_Research_Report_EN.pdf",
            "FIN_Local_Research_Checkpoint_P428_P430_State.json",
            "FIN_Programs_428_430_Results.json",
            "RELEASE_10_36_PROGRAMS_428_430.md",
            "FIN_Checkpoint_P428_P430_AGENTS_Guardrail.txt",
        ):
            self.assertTrue((ROOT / name).is_file(), name)

    def test_checkpoint_state_preserves_boundaries(self) -> None:
        state = json.loads(
            (ROOT / "FIN_Local_Research_Checkpoint_P428_P430_State.json").read_text(
                encoding="utf-8"
            )
        )
        self.assertEqual(state["latest_completed_program"], "P430")
        self.assertEqual(state["next_selected_batch"], ["P435", "P436", "P440"])
        self.assertTrue(state["p429_exact_krawczyk_inclusion"])
        self.assertTrue(state["p430_global_dual_feasibility"])
        self.assertFalse(state["physical_evidence_admitted"])
        self.assertFalse(state["selector_discharged"])
        self.assertFalse(state["dimensional_source_exported"])
        self.assertFalse(state["legacy_strict_bridge_complete"])
        self.assertFalse(state["legacy_star_typed_self_coupling_exported"])

    def test_report_contains_required_scientific_boundaries(self) -> None:
        text = (ROOT / "FIN_Local_Research_Checkpoint_P428_P430_EN.md").read_text(
            encoding="utf-8"
        )
        for token in (
            "QW-2191",
            "Rational-Cut Cosine Provider Interface",
            "Exact-Rational Seven-Contact KKT Isolating Box",
            "Contact-Aware Global Dual Certificate",
            "K^*_{\\mathrm{legacy}}",
            "not yet a mathematical self-coupling law",
            "No complete legacy-to-strict bridge",
            "Restart instructions",
        ):
            self.assertIn(token, text)

    def test_pdf_metadata_and_text(self) -> None:
        pdf = ROOT / "FIN_Local_Research_Checkpoint_P428_P430_EN.pdf"
        self.assertGreater(pdf.stat().st_size, 400_000)
        info = subprocess.run(
            ["pdfinfo", pdf.name],
            cwd=ROOT,
            capture_output=True,
            text=True,
            timeout=30,
            check=True,
        ).stdout
        self.assertIn("FIN Local Research Checkpoint P428", info)
        self.assertIn("Krzysztof Żuchowski", info)
        pages = next(
            int(line.split(":", 1)[1].strip())
            for line in info.splitlines()
            if line.startswith("Pages:")
        )
        self.assertGreaterEqual(pages, 12)

        extracted = subprocess.run(
            ["pdftotext", "-layout", pdf.name, "-"],
            cwd=ROOT,
            capture_output=True,
            text=True,
            timeout=30,
            check=True,
        ).stdout
        for token in (
            "Krzysztof Żuchowski",
            "4.0672647472032875",
            "0.7073534677231137",
            "P435/P436/P440",
            "No laboratory data",
        ):
            self.assertIn(token, extracted)

    def test_checkpoint_is_retained_after_rolling_report_advances(self) -> None:
        checkpoint = ROOT / "FIN_Local_Research_Checkpoint_P428_P430_EN.pdf"
        current = ROOT / "FIN_Current_Local_Research_Report_EN.pdf"
        self.assertTrue(checkpoint.is_file())
        self.assertGreater(checkpoint.stat().st_size, 400_000)
        candidates = (
            ROOT / "FIN_Local_Research_Checkpoint_P458_P464_EN.pdf",
            ROOT / "FIN_Local_Research_Checkpoint_P454_P457_EN.pdf",
            ROOT / "FIN_Local_Research_Checkpoint_P451_P453_EN.pdf",
            ROOT / "FIN_Local_Research_Checkpoint_P448_P450_EN.pdf",
            ROOT / "FIN_Local_Research_Checkpoint_P445_P447_EN.pdf",
            ROOT / "FIN_Local_Research_Checkpoint_P435_P440_EN.pdf",
            checkpoint,
        )
        latest = next(path for path in candidates if path.is_file())
        self.assertEqual(latest.read_bytes(), current.read_bytes())


if __name__ == "__main__":
    unittest.main()
