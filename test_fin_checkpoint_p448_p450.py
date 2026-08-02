#!/usr/bin/env python3
"""Integrity tests for the FIN Release-10.39 durable checkpoint."""

from __future__ import annotations

import json
from pathlib import Path
import subprocess
import unittest


ROOT = Path(__file__).resolve().parent


class CheckpointP448P450Tests(unittest.TestCase):
    def test_required_artifacts_exist(self) -> None:
        for name in (
            "FIN_Local_Research_Checkpoint_P448_P450_EN.md",
            "FIN_Local_Research_Checkpoint_P448_P450_EN.tex",
            "FIN_Local_Research_Checkpoint_P448_P450_EN.pdf",
            "FIN_Current_Local_Research_Report_EN.pdf",
            "FIN_Local_Research_Checkpoint_P448_P450_State.json",
            "FIN_Programs_448_450_Results.json",
            "FIN_Programs_448_450_P449_Three_Slot_Witness.npz",
            "RELEASE_10_39_PROGRAMS_448_450.md",
            "FIN_Checkpoint_P448_P450_AGENTS_Guardrail.txt",
        ):
            self.assertTrue((ROOT / name).is_file(), name)

    def test_state_preserves_results_and_boundaries(self) -> None:
        state = json.loads(
            (ROOT / "FIN_Local_Research_Checkpoint_P448_P450_State.json").read_text(
                encoding="utf-8"
            )
        )
        self.assertEqual(state["latest_completed_program"], "P450")
        self.assertEqual(state["next_selected_batch"], ["P451", "P452", "P453"])
        self.assertTrue(state["p448_full_simplex_global_upper_proved"])
        self.assertFalse(state["p448_exact_full_simplex_globality_proved"])
        self.assertTrue(state["p449_recursive_support_cone_proved"])
        self.assertTrue(state["p449_echo_history_full_extension_exported"])
        self.assertTrue(state["p449_echo_beats_ghz_certified"])
        self.assertGreater(state["p449_echo_advantage_lower"], 0.0176)
        self.assertFalse(state["p449_global_three_slot_optimum_proved"])
        self.assertTrue(state["p450_representation_independence_refuted"])
        self.assertFalse(state["p450_canonical_representation_exported"])
        self.assertFalse(state["physical_evidence_admitted"])
        self.assertFalse(state["selector_discharged"])
        self.assertFalse(state["dimensional_source_exported"])
        self.assertFalse(state["legacy_strict_bridge_complete"])

    def test_report_contains_required_results_and_nonclaims(self) -> None:
        text = (ROOT / "FIN_Local_Research_Checkpoint_P448_P450_EN.md").read_text(
            encoding="utf-8"
        )
        for token in (
            "O164 -- Concave Fine-Grained-Erasure Majorant",
            "O165 -- Recursive Causal-Support Cone",
            "O166 -- Moment-Null Representation Gauge",
            "0.4666305033804779",
            "0.49063899018433244",
            "0.017612669538553616",
            "Three-slot GHZ optimality is therefore **refuted**",
            "sup_\\epsilon V(\\mu_\\epsilon)=+\\infty",
            "QW-2191",
            "P451\\longrightarrow P452\\longrightarrow P453",
            "restart instructions",
        ):
            self.assertIn(token, text)

    def test_pdf_metadata_text_figure_and_clean_log(self) -> None:
        pdf = ROOT / "FIN_Local_Research_Checkpoint_P448_P450_EN.pdf"
        self.assertGreater(pdf.stat().st_size, 480_000)
        info = subprocess.run(
            ["pdfinfo", pdf.name], cwd=ROOT, capture_output=True, text=True,
            timeout=30, check=True,
        ).stdout
        self.assertIn("FIN Local Research Checkpoint P448", info)
        self.assertIn("Krzysztof Żuchowski", info)
        pages = next(
            int(line.split(":", 1)[1].strip())
            for line in info.splitlines() if line.startswith("Pages:")
        )
        self.assertGreaterEqual(pages, 28)

        extracted = subprocess.run(
            ["pdftotext", "-layout", pdf.name, "-"], cwd=ROOT,
            capture_output=True, text=True, timeout=30, check=True,
        ).stdout
        for token in (
            "Krzysztof Żuchowski",
            "0.4666305033804779",
            "0.49063899018433244",
            "0.017612669538553616",
            "367.38655583",
            "P451",
        ):
            self.assertIn(token, extracted)
        self.assertGreater(
            (ROOT / "FIN_Programs_448_450_Figures/p448_p450_global_majorant_and_gauge_obstruction.png").stat().st_size,
            120_000,
        )
        log = (ROOT / "FIN_Local_Research_Checkpoint_P448_P450_EN.log").read_text(
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
        )
        checkpoint_path = next(path for path in candidates if path.is_file())
        checkpoint = checkpoint_path.read_bytes()
        current = (ROOT / "FIN_Current_Local_Research_Report_EN.pdf").read_bytes()
        self.assertEqual(checkpoint, current)


if __name__ == "__main__":
    unittest.main()
