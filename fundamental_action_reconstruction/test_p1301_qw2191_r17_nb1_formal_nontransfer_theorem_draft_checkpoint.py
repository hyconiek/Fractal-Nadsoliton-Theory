from __future__ import annotations

import json
import subprocess
import tempfile
import unittest
from pathlib import Path

SCRIPT = Path(__file__).resolve().parent / "p1301_qw2191_r17_nb1_formal_nontransfer_theorem_draft_checkpoint.py"


class TestP1301QW2191R17NB1FormalNontransferTheoremDraftCheckpoint(unittest.TestCase):
    def test_happy_path(self) -> None:
        with tempfile.TemporaryDirectory() as td:
            td_path = Path(td)
            p1300 = td_path / "p1300.json"
            out = td_path / "p1301.json"
            p1300.write_text(
                json.dumps(
                    {
                        "next_priority": "R17_NB1_FORMAL_NONTRANSFER_THEOREM_DRAFT",
                        "r16_nonbridge_scope_limit_theorem": {"status": "SCOPE_LIMIT_DRAFTED"},
                    }
                ),
                encoding="utf-8",
            )
            subprocess.run(["python3", str(SCRIPT), "--p1300", str(p1300), "--out", str(out)], check=True)
            payload = json.loads(out.read_text(encoding="utf-8"))
            self.assertEqual(payload["r17_nb1_nontransfer_theorem"]["status"], "DRAFT_WITH_OBLIGATIONS")

    def test_requires_r16_status(self) -> None:
        with tempfile.TemporaryDirectory() as td:
            td_path = Path(td)
            p1300 = td_path / "p1300.json"
            out = td_path / "p1301.json"
            p1300.write_text(
                json.dumps(
                    {
                        "next_priority": "R17_NB1_FORMAL_NONTRANSFER_THEOREM_DRAFT",
                        "r16_nonbridge_scope_limit_theorem": {"status": "BLOCKED"},
                    }
                ),
                encoding="utf-8",
            )
            proc = subprocess.run(["python3", str(SCRIPT), "--p1300", str(p1300), "--out", str(out)], capture_output=True, text=True)
            self.assertNotEqual(proc.returncode, 0)
            self.assertIn("SCOPE_LIMIT_DRAFTED", proc.stderr)


if __name__ == "__main__":
    unittest.main()
