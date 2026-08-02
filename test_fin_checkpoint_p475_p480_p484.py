#!/usr/bin/env python3
"""Integrity tests for the FIN Release-10.46 checkpoint."""

import json
from pathlib import Path
import subprocess
import unittest


ROOT = Path(__file__).resolve().parent


class CheckpointP475P480P484Tests(unittest.TestCase):
    def test_required_artifacts(self) -> None:
        for name in (
            "FIN_Local_Research_Checkpoint_P475_P480_P484_EN.md",
            "FIN_Local_Research_Checkpoint_P475_P480_P484_EN.tex",
            "FIN_Local_Research_Checkpoint_P475_P480_P484_EN.pdf",
            "FIN_Local_Research_Checkpoint_P475_P480_P484_State.json",
            "FIN_Program_475_Results.json",
            "FIN_Program_475_Algebraic_Inventory.csv",
            "FIN_Program_480_Standalone_Certificate.json",
            "FIN_Program_480_Standalone_Check_Result.json",
            "FIN_Program_484_Results.json",
            "FIN_Program_484_Phase_Face_Witness.npz",
            "RELEASE_10_46_PROGRAMS_475_480_484.md",
            "FIN_Checkpoint_P475_P480_P484_AGENTS_Guardrail.txt",
            "FIN_PROGRAMS_475_480_484_RELEASE_MANIFEST.sha256",
        ):
            self.assertTrue((ROOT / name).is_file(), name)

    def test_state_and_boundaries(self) -> None:
        state = json.loads(
            (ROOT / "FIN_Local_Research_Checkpoint_P475_P480_P484_State.json").read_text(
                encoding="utf-8"
            )
        )
        self.assertEqual(state["next_batch"], ["P485", "P486", "P487"])
        self.assertEqual(state["programs"]["P475"]["status"],
                         "open_resource_bounded_no_go")
        self.assertEqual(state["programs"]["P480"]["status"],
                         "computer_assisted_proof_replayed")
        self.assertEqual(state["programs"]["P484"]["status"],
                         "strong_evidence_exact_causal_axis_identity_open")
        self.assertEqual(state["selector_QW_2191"], "open")
        self.assertEqual(state["dimensional_source"], "missing")
        self.assertFalse(state["network_used"])
        self.assertFalse(state["firecrawl_used"])
        self.assertFalse(state["laboratory_data_used"])
        self.assertFalse(state["external_audit_used"])

    def test_pdf(self) -> None:
        pdf = ROOT / "FIN_Local_Research_Checkpoint_P475_P480_P484_EN.pdf"
        self.assertGreater(pdf.stat().st_size, 400_000)
        info = subprocess.run(
            ["pdfinfo", pdf.name], cwd=ROOT, capture_output=True, text=True, check=True
        ).stdout
        self.assertIn("Krzysztof Żuchowski", info)
        self.assertIn("A4", info)
        self.assertIn("Pages:           19", info)
        extracted = subprocess.run(
            ["pdftotext", "-layout", pdf.name, "-"],
            cwd=ROOT,
            capture_output=True,
            text=True,
            check=True,
        ).stdout
        for token in (
            "P475",
            "P480",
            "P484",
            "O185",
            "O186",
            "O187",
            "QW-2191",
            "P485",
            "resource-bounded no-go",
            "Full-cone optimizer uniqueness is neither proved nor refuted",
        ):
            self.assertIn(token, extracted)
        log = (ROOT / "FIN_Local_Research_Checkpoint_P475_P480_P484_EN.log").read_text(
            encoding="utf-8", errors="replace"
        )
        self.assertNotIn("LaTeX Error", log)
        self.assertNotIn("Missing character", log)
        self.assertNotIn("Fatal error", log)
        self.assertNotIn("Overfull", log)

    def test_current_copy(self) -> None:
        self.assertEqual(
            (ROOT / "FIN_Local_Research_Checkpoint_P475_P480_P484_EN.pdf").read_bytes(),
            (ROOT / "FIN_Current_Local_Research_Report_EN.pdf").read_bytes(),
        )

    def test_release_manifest(self) -> None:
        completed = subprocess.run(
            ["sha256sum", "-c", "FIN_PROGRAMS_475_480_484_RELEASE_MANIFEST.sha256"],
            cwd=ROOT,
            capture_output=True,
            text=True,
            check=True,
        )
        self.assertEqual(completed.stdout.count(": OK"), 24)


if __name__ == "__main__":
    unittest.main()
