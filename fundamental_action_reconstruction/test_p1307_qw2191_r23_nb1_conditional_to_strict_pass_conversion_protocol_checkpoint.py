from __future__ import annotations

import json
import subprocess
import tempfile
import unittest
from pathlib import Path

SCRIPT = Path(__file__).resolve().parent / "p1307_qw2191_r23_nb1_conditional_to_strict_pass_conversion_protocol_checkpoint.py"


class TestP1307QW2191R23NB1ConditionalToStrictPassConversionProtocolCheckpoint(unittest.TestCase):
    def test_happy_path(self) -> None:
        with tempfile.TemporaryDirectory() as td:
            td_path = Path(td)
            p1306 = td_path / "p1306.json"
            out = td_path / "p1307.json"
            p1306.write_text(
                json.dumps(
                    {
                        "next_priority": "R23_NB1_CONDITIONAL_TO_STRICT_PASS_CONVERSION_PROTOCOL",
                        "batch_status": "NEAR_COMPLETE",
                    }
                ),
                encoding="utf-8",
            )
            subprocess.run(["python3", str(SCRIPT), "--p1306", str(p1306), "--out", str(out)], check=True)
            payload = json.loads(out.read_text(encoding="utf-8"))
            self.assertEqual(payload["r23_conversion_protocol"]["status"], "PROTOCOL_READY")

    def test_requires_near_complete(self) -> None:
        with tempfile.TemporaryDirectory() as td:
            td_path = Path(td)
            p1306 = td_path / "p1306.json"
            out = td_path / "p1307.json"
            p1306.write_text(
                json.dumps(
                    {
                        "next_priority": "R23_NB1_CONDITIONAL_TO_STRICT_PASS_CONVERSION_PROTOCOL",
                        "batch_status": "PARTIAL_PROGRESS",
                    }
                ),
                encoding="utf-8",
            )
            proc = subprocess.run(["python3", str(SCRIPT), "--p1306", str(p1306), "--out", str(out)], capture_output=True, text=True)
            self.assertNotEqual(proc.returncode, 0)
            self.assertIn("NEAR_COMPLETE", proc.stderr)


if __name__ == "__main__":
    unittest.main()
