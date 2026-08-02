#!/usr/bin/env python3
"""Integrity tests for FIN goal-mode Checkpoint 000."""

from __future__ import annotations

import json
from pathlib import Path
import unittest


ROOT = Path(__file__).resolve().parent


class GoalCheckpoint000Tests(unittest.TestCase):
    def test_required_sources_exist(self) -> None:
        for name in (
            "FIN_Local_Research_Checkpoint_000_EN.md",
            "FIN_Local_Research_Checkpoint_000_EN.tex",
            "FIN_Local_Research_Checkpoint_000_State.json",
            "FIN_GOAL_STATE.md",
            "RELEASE_10_36_CHECKPOINT_000.md",
        ):
            self.assertTrue((ROOT / name).is_file(), name)

    def test_state_preserves_boundaries(self) -> None:
        state = json.loads(
            (ROOT / "FIN_Local_Research_Checkpoint_000_State.json").read_text(
                encoding="utf-8"
            )
        )
        self.assertEqual(state["latest_completed_program"], "P427")
        self.assertEqual(state["selected_batch"], ["P428", "P429", "P430"])
        self.assertFalse(state["physical_evidence_admitted"])
        self.assertFalse(state["selector_discharged"])
        self.assertFalse(state["dimensional_source_exported"])
        self.assertFalse(state["legacy_strict_bridge_complete"])
        self.assertFalse(state["role_transfer_started"])

    def test_report_contains_required_boundaries(self) -> None:
        text = (ROOT / "FIN_Local_Research_Checkpoint_000_EN.md").read_text(
            encoding="utf-8"
        )
        for token in (
            "QW-2191",
            "P428",
            "P429",
            "P430",
            "Synthetic records",
            "Explicit nonclaims",
            "Restart instructions",
            "legacy-to-strict",
        ):
            self.assertIn(token, text)

    def test_pdf_is_durable_when_built(self) -> None:
        pdf = ROOT / "FIN_Local_Research_Checkpoint_000_EN.pdf"
        if pdf.exists():
            self.assertGreater(pdf.stat().st_size, 50_000)


if __name__ == "__main__":
    unittest.main()

