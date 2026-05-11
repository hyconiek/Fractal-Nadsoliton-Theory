from __future__ import annotations

import json
import subprocess
import tempfile
import unittest
from pathlib import Path

SCRIPT = Path(__file__).resolve().parent / "p1292_qw2191_r8b_selector_split_rerun_with_updated_margin_checkpoint.py"


class TestP1292QW2191R8BSelectorSplitRerunWithUpdatedMarginCheckpoint(unittest.TestCase):
    def test_happy_path(self) -> None:
        with tempfile.TemporaryDirectory() as td:
            td_path = Path(td)
            p1291 = td_path / "p1291.json"
            out = td_path / "p1292.json"
            p1291.write_text(
                json.dumps(
                    {
                        "next_priority": "R8B_SELECTOR_SPLIT_RERUN_WITH_UPDATED_MARGIN",
                        "r8a_execution": {"acceptance_pass": True},
                    }
                ),
                encoding="utf-8",
            )
            subprocess.run(["python3", str(SCRIPT), "--p1291", str(p1291), "--out", str(out)], check=True)
            payload = json.loads(out.read_text(encoding="utf-8"))
            self.assertEqual(payload["r8b_split_rerun"]["report"]["result"], "DECISIVE_FOR_SSEL_SRC_A")
            self.assertEqual(payload["next_priority"], "R9_FORMAL_SELECTOR_SOURCE_THEOREM_DRAFT")

    def test_requires_passed_r8a(self) -> None:
        with tempfile.TemporaryDirectory() as td:
            td_path = Path(td)
            p1291 = td_path / "p1291.json"
            out = td_path / "p1292.json"
            p1291.write_text(
                json.dumps(
                    {
                        "next_priority": "R8B_SELECTOR_SPLIT_RERUN_WITH_UPDATED_MARGIN",
                        "r8a_execution": {"acceptance_pass": False},
                    }
                ),
                encoding="utf-8",
            )
            proc = subprocess.run(["python3", str(SCRIPT), "--p1291", str(p1291), "--out", str(out)], capture_output=True, text=True)
            self.assertNotEqual(proc.returncode, 0)
            self.assertIn("acceptance_pass=true", proc.stderr)


if __name__ == "__main__":
    unittest.main()
