from __future__ import annotations

import json
import subprocess
import tempfile
import unittest
from pathlib import Path

SCRIPT = Path(__file__).resolve().parent / "p1309_qw2191_r25_nb1_formal_closure_statement_and_export_packet_checkpoint.py"


class TestP1309QW2191R25NB1FormalClosureStatementAndExportPacketCheckpoint(unittest.TestCase):
    def test_happy_path(self) -> None:
        with tempfile.TemporaryDirectory() as td:
            td_path = Path(td)
            p1308 = td_path / "p1308.json"
            out = td_path / "p1309.json"
            p1308.write_text(
                json.dumps(
                    {
                        "next_priority": "R25_NB1_FORMAL_CLOSURE_STATEMENT_AND_EXPORT_PACKET",
                        "lnb1_2_final_status": "PASS_STRICT",
                    }
                ),
                encoding="utf-8",
            )
            subprocess.run(["python3", str(SCRIPT), "--p1308", str(p1308), "--out", str(out)], check=True)
            payload = json.loads(out.read_text(encoding="utf-8"))
            self.assertEqual(payload["r25_closure_statement"]["status"], "FORMAL_CLOSURE_DECLARED")

    def test_requires_pass_strict(self) -> None:
        with tempfile.TemporaryDirectory() as td:
            td_path = Path(td)
            p1308 = td_path / "p1308.json"
            out = td_path / "p1309.json"
            p1308.write_text(
                json.dumps(
                    {
                        "next_priority": "R25_NB1_FORMAL_CLOSURE_STATEMENT_AND_EXPORT_PACKET",
                        "lnb1_2_final_status": "PASS_CONDITIONAL",
                    }
                ),
                encoding="utf-8",
            )
            proc = subprocess.run(["python3", str(SCRIPT), "--p1308", str(p1308), "--out", str(out)], capture_output=True, text=True)
            self.assertNotEqual(proc.returncode, 0)
            self.assertIn("PASS_STRICT", proc.stderr)


if __name__ == "__main__":
    unittest.main()
