#!/usr/bin/env python3
"""Integrity tests for the FIN Release-10.44 checkpoint."""

import json
from pathlib import Path
import subprocess
import unittest


ROOT = Path(__file__).resolve().parent


class CheckpointP471P473Tests(unittest.TestCase):
    def test_required_artifacts(self) -> None:
        for name in (
            "FIN_Local_Research_Checkpoint_P471_P473_EN.md",
            "FIN_Local_Research_Checkpoint_P471_P473_EN.tex",
            "FIN_Local_Research_Checkpoint_P471_P473_EN.pdf",
            "FIN_Local_Research_Checkpoint_P471_P473_State.json",
            "FIN_Programs_471_472_473_Results.json",
            "FIN_Programs_471_472_473_P473_Root_Box.npz",
            "RELEASE_10_44_PROGRAMS_471_472_473.md",
            "FIN_Checkpoint_P471_P473_AGENTS_Guardrail.txt",
        ):
            self.assertTrue((ROOT / name).is_file(), name)

    def test_state_and_boundaries(self) -> None:
        state = json.loads(
            (ROOT / "FIN_Local_Research_Checkpoint_P471_P473_State.json").read_text(
                encoding="utf-8"
            )
        )
        self.assertEqual(state["next_selected_batch"], ["P474", "P483", "P479"])
        self.assertTrue(state["p472_strict_krawczyk_inclusion"])
        self.assertTrue(state["p473_exact_root_proved"])
        self.assertTrue(state["p473_exact_O167_attainment_proved"])
        self.assertFalse(state["p473_full_cone_optimizer_uniqueness_proved"])
        self.assertLess(state["p473_global_value_width"], 1e-13)
        for key in (
            "physical_evidence_admitted",
            "selector_discharged",
            "dimensional_source_exported",
            "legacy_strict_bridge_complete",
            "legacy_role_transfer_started",
        ):
            self.assertFalse(state[key])

    def test_pdf(self) -> None:
        pdf = ROOT / "FIN_Local_Research_Checkpoint_P471_P473_EN.pdf"
        self.assertGreater(pdf.stat().st_size, 420_000)
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
            "0.5233281002710717",
            "O179",
            "O180",
            "O181",
            "P474",
            "globally optimal",
        ):
            self.assertIn(token, extracted)
        log = (ROOT / "FIN_Local_Research_Checkpoint_P471_P473_EN.log").read_text(
            encoding="utf-8", errors="replace"
        )
        self.assertNotIn("LaTeX Error", log)
        self.assertNotIn("Missing character", log)
        self.assertNotIn("Fatal error", log)
        self.assertNotIn("Overfull", log)

    def test_current_copy(self) -> None:
        self.assertEqual(
            (ROOT / "FIN_Local_Research_Checkpoint_P471_P473_EN.pdf").read_bytes(),
            (ROOT / "FIN_Current_Local_Research_Report_EN.pdf").read_bytes(),
        )


if __name__ == "__main__":
    unittest.main()
