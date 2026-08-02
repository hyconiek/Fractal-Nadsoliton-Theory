#!/usr/bin/env python3
"""Integrity tests for the FIN Release-10.41 durable checkpoint."""

from __future__ import annotations

import json
from pathlib import Path
import subprocess
import unittest


ROOT = Path(__file__).resolve().parent


class CheckpointP454P457Tests(unittest.TestCase):
    def test_required_artifacts_exist(self) -> None:
        for name in (
            "FIN_Local_Research_Checkpoint_P454_P457_EN.md",
            "FIN_Local_Research_Checkpoint_P454_P457_EN.tex",
            "FIN_Local_Research_Checkpoint_P454_P457_EN.pdf",
            "FIN_Current_Local_Research_Report_EN.pdf",
            "FIN_Local_Research_Checkpoint_P454_P457_State.json",
            "FIN_Programs_454_455_457_Results.json",
            "FIN_Programs_454_455_457_P454_Rational_Dual.npz",
            "FIN_Programs_454_455_457_P454_Nested_Dual.csv",
            "FIN_Programs_454_455_457_P455_Symmetry_Residual.csv",
            "FIN_Programs_454_455_457_P457_Refined_Cover.csv",
            "RELEASE_10_41_PROGRAMS_454_455_457.md",
            "FIN_Checkpoint_P454_P457_AGENTS_Guardrail.txt",
        ):
            self.assertTrue((ROOT / name).is_file(), name)

    def test_state_preserves_results_and_boundaries(self) -> None:
        state = json.loads(
            (ROOT / "FIN_Local_Research_Checkpoint_P454_P457_State.json").read_text(
                encoding="utf-8"
            )
        )
        self.assertEqual(state["release"], "10.41")
        self.assertEqual(state["latest_completed_program"], "P457")
        self.assertEqual(state["next_selected_batch"], ["P464", "P458", "P459"])
        self.assertTrue(state["p454_nested_dual_derived"])
        self.assertTrue(state["p454_rational_dual_feasible"])
        self.assertLess(state["p454_global_gap"], 6.7e-6)
        self.assertFalse(state["p454_exact_optimizer_proved"])
        self.assertFalse(state["p454_optimizer_unique"])
        self.assertEqual(state["p455_real_complement_fixed_dimension"], 5)
        self.assertEqual(state["p455_o167_face_dimension"], 3)
        self.assertTrue(state["p455_symmetry_shortcut_refuted"])
        self.assertFalse(state["p455_symmetry_reduction_to_o167_proved"])
        self.assertLess(state["p457_full_simplex_gap"], 1e-6)
        self.assertFalse(state["physical_evidence_admitted"])
        self.assertFalse(state["selector_discharged"])
        self.assertFalse(state["dimensional_source_exported"])
        self.assertFalse(state["legacy_strict_bridge_complete"])
        self.assertFalse(state["legacy_role_transfer_started"])

    def test_report_contains_required_results_and_nonclaims(self) -> None:
        text = (ROOT / "FIN_Local_Research_Checkpoint_P454_P457_EN.md").read_text(
            encoding="utf-8"
        )
        for token in (
            "O170 — Nested Comb Dual Ladder",
            "O171 — Ordered-Comb Symmetry Residual Space",
            "O172 — Six-Decimal Coarse-Erasure Cover",
            "0.52332810026048937",
            "0.523334700252",
            "4.9999999969",
            "five-dimensional fixed space",
            "two symmetry-allowed directions",
            "0.46327928294340853",
            "QW-2191 | Open",
            "P464\\longrightarrow P458\\longrightarrow P459",
            "Restart instructions",
        ):
            self.assertIn(token, text)

    def test_pdf_metadata_text_figure_and_clean_log(self) -> None:
        pdf = ROOT / "FIN_Local_Research_Checkpoint_P454_P457_EN.pdf"
        self.assertGreater(pdf.stat().st_size, 490_000)
        info = subprocess.run(
            ["pdfinfo", pdf.name], cwd=ROOT, capture_output=True, text=True,
            timeout=30, check=True,
        ).stdout
        self.assertIn("FIN Local Research Checkpoint P454", info)
        self.assertIn("Krzysztof Żuchowski", info)
        self.assertIn("A4", info)
        pages = next(
            int(line.split(":", 1)[1].strip())
            for line in info.splitlines() if line.startswith("Pages:")
        )
        self.assertGreaterEqual(pages, 27)

        extracted = subprocess.run(
            ["pdftotext", "-layout", pdf.name, "-"], cwd=ROOT,
            capture_output=True, text=True, timeout=30, check=True,
        ).stdout
        for token in (
            "Krzysztof Żuchowski",
            "0.52332810026048937",
            "0.523334700252",
            "O170",
            "O171",
            "O172",
            "P464",
        ):
            self.assertIn(token, extracted)
        self.assertGreater(
            (ROOT / "FIN_Programs_454_455_457_Figures/p454_p457_dual_symmetry_and_refined_cover.png").stat().st_size,
            160_000,
        )
        log = (ROOT / "FIN_Local_Research_Checkpoint_P454_P457_EN.log").read_text(
            encoding="utf-8", errors="replace"
        )
        self.assertNotIn("Missing character", log)
        self.assertNotIn("LaTeX Error", log)
        self.assertNotIn("Undefined control sequence", log)

    def test_current_report_matches_latest_checkpoint(self) -> None:
        latest = ROOT / "FIN_Local_Research_Checkpoint_P458_P464_EN.pdf"
        checkpoint_path = latest if latest.is_file() else ROOT / "FIN_Local_Research_Checkpoint_P454_P457_EN.pdf"
        checkpoint = checkpoint_path.read_bytes()
        current = (ROOT / "FIN_Current_Local_Research_Report_EN.pdf").read_bytes()
        self.assertEqual(checkpoint, current)


if __name__ == "__main__":
    unittest.main()
