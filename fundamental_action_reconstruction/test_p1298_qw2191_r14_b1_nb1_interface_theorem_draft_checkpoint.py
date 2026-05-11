from __future__ import annotations

import json
import subprocess
import tempfile
import unittest
from pathlib import Path

SCRIPT = Path(__file__).resolve().parent / "p1298_qw2191_r14_b1_nb1_interface_theorem_draft_checkpoint.py"


class TestP1298QW2191R14B1NB1InterfaceTheoremDraftCheckpoint(unittest.TestCase):
    def test_happy_path(self) -> None:
        with tempfile.TemporaryDirectory() as td:
            td_path = Path(td)
            p1297 = td_path / "p1297.json"
            out = td_path / "p1298.json"
            p1297.write_text(
                json.dumps(
                    {
                        "next_priority": "R14_B1_NB1_INTERFACE_THEOREM_DRAFT",
                        "r13_interface": {"status": "DRAFT_READY"},
                    }
                ),
                encoding="utf-8",
            )
            subprocess.run(["python3", str(SCRIPT), "--p1297", str(p1297), "--out", str(out)], check=True)
            payload = json.loads(out.read_text(encoding="utf-8"))
            self.assertEqual(payload["r14_interface_theorem"]["status"], "THEOREM_INTERFACE_DRAFTED")
            self.assertEqual(payload["next_priority"], "R15_B1_NB1_OBLIGATION_MATRIX_AND_PROOF_PLAN")

    def test_requires_draft_ready(self) -> None:
        with tempfile.TemporaryDirectory() as td:
            td_path = Path(td)
            p1297 = td_path / "p1297.json"
            out = td_path / "p1298.json"
            p1297.write_text(
                json.dumps(
                    {
                        "next_priority": "R14_B1_NB1_INTERFACE_THEOREM_DRAFT",
                        "r13_interface": {"status": "BLOCKED"},
                    }
                ),
                encoding="utf-8",
            )
            proc = subprocess.run(["python3", str(SCRIPT), "--p1297", str(p1297), "--out", str(out)], capture_output=True, text=True)
            self.assertNotEqual(proc.returncode, 0)
            self.assertIn("DRAFT_READY", proc.stderr)


if __name__ == "__main__":
    unittest.main()
