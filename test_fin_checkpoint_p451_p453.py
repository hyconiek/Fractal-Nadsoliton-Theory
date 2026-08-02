#!/usr/bin/env python3
"""Integrity tests for the FIN Release-10.40 durable checkpoint."""

from __future__ import annotations

import json
from pathlib import Path
import subprocess
import unittest


ROOT = Path(__file__).resolve().parent


class CheckpointP451P453Tests(unittest.TestCase):
    def test_required_artifacts_exist(self) -> None:
        for name in (
            "FIN_Local_Research_Checkpoint_P451_P453_EN.md",
            "FIN_Local_Research_Checkpoint_P451_P453_EN.tex",
            "FIN_Local_Research_Checkpoint_P451_P453_EN.pdf",
            "FIN_Current_Local_Research_Report_EN.pdf",
            "FIN_Local_Research_Checkpoint_P451_P453_State.json",
            "FIN_Programs_451_453_Results.json",
            "FIN_Programs_451_453_P451_Coherent_Witness.npz",
            "RELEASE_10_40_PROGRAMS_451_453.md",
            "FIN_Checkpoint_P451_P453_AGENTS_Guardrail.txt",
        ):
            self.assertTrue((ROOT / name).is_file(), name)

    def test_state_preserves_results_and_boundaries(self) -> None:
        state = json.loads(
            (ROOT / "FIN_Local_Research_Checkpoint_P451_P453_State.json").read_text(
                encoding="utf-8"
            )
        )
        self.assertEqual(state["latest_completed_program"], "P453")
        self.assertEqual(state["next_selected_batch"], ["P454", "P455", "P457"])
        self.assertTrue(state["p451_global_diagonal_optimum_proved"])
        self.assertTrue(state["p451_coherent_advantage_certified"])
        self.assertGreater(state["p451_coherent_advantage_lower"], 5.6e-4)
        self.assertFalse(state["p451_global_full_cone_optimum_proved"])
        self.assertTrue(state["p452_coarse_reversal_symmetrization_proved"])
        self.assertTrue(state["p452_full_simplex_certificate_proved"])
        self.assertLess(state["p452_full_simplex_gap"], 1e-3)
        self.assertTrue(state["p453_global_minimum_negativity_representation_unique"])
        self.assertFalse(state["p453_moments_alone_canonical"])
        self.assertFalse(state["physical_evidence_admitted"])
        self.assertFalse(state["selector_discharged"])
        self.assertFalse(state["dimensional_source_exported"])
        self.assertFalse(state["legacy_strict_bridge_complete"])

    def test_report_contains_required_results_and_nonclaims(self) -> None:
        text = (ROOT / "FIN_Local_Research_Checkpoint_P451_P453_EN.md").read_text(
            encoding="utf-8"
        )
        for token in (
            "O167 — Causal-Coherence Advantage Witness",
            "O168 — Coarse-Erasure Symmetrization Certificate",
            "O169 — Strict-Complementarity Jordan Gauge Fixer",
            "0.52332810026048937",
            "0.00056473095851463685",
            "0.4642707484333634",
            "1.4954513576571674",
            "diagonal optimality is **[Refuted]**",
            "QW-2191 remains open",
            "P454\\longrightarrow P455\\longrightarrow P457",
            "Restart instructions",
        ):
            self.assertIn(token, text)

    def test_pdf_metadata_text_figure_and_clean_log(self) -> None:
        pdf = ROOT / "FIN_Local_Research_Checkpoint_P451_P453_EN.pdf"
        self.assertGreater(pdf.stat().st_size, 530_000)
        info = subprocess.run(
            ["pdfinfo", pdf.name], cwd=ROOT, capture_output=True, text=True,
            timeout=30, check=True,
        ).stdout
        self.assertIn("FIN Local Research Checkpoint P451", info)
        self.assertIn("Krzysztof Żuchowski", info)
        self.assertIn("A4", info)
        pages = next(
            int(line.split(":", 1)[1].strip())
            for line in info.splitlines() if line.startswith("Pages:")
        )
        self.assertGreaterEqual(pages, 31)

        extracted = subprocess.run(
            ["pdftotext", "-layout", pdf.name, "-"], cwd=ROOT,
            capture_output=True, text=True, timeout=30, check=True,
        ).stdout
        for token in (
            "Krzysztof Żuchowski",
            "0.52332810026048937",
            "0.00056473095851463685",
            "0.4642707484333634",
            "1.4954513576571674",
            "P454",
        ):
            self.assertIn(token, extracted)
        self.assertGreater(
            (ROOT / "FIN_Programs_451_453_Figures/p451_p453_coherence_symmetry_and_gauge_fixing.png").stat().st_size,
            170_000,
        )
        log = (ROOT / "FIN_Local_Research_Checkpoint_P451_P453_EN.log").read_text(
            encoding="utf-8", errors="replace"
        )
        self.assertNotIn("Missing character", log)
        self.assertNotIn("LaTeX Error", log)
        self.assertNotIn("Undefined control sequence", log)

    def test_current_report_matches_latest_checkpoint(self) -> None:
        candidates = (
            ROOT / "FIN_Local_Research_Checkpoint_P458_P464_EN.pdf",
            ROOT / "FIN_Local_Research_Checkpoint_P454_P457_EN.pdf",
            ROOT / "FIN_Local_Research_Checkpoint_P451_P453_EN.pdf",
        )
        checkpoint_path = next(path for path in candidates if path.is_file())
        checkpoint = checkpoint_path.read_bytes()
        current = (ROOT / "FIN_Current_Local_Research_Report_EN.pdf").read_bytes()
        self.assertEqual(checkpoint, current)


if __name__ == "__main__":
    unittest.main()
