from __future__ import annotations

import json
import subprocess
import tempfile
import unittest
from pathlib import Path

SCRIPT = Path(__file__).resolve().parent / "p1290_qw2191_r8a_data_acquisition_protocol_v2_checkpoint.py"


class TestP1290QW2191R8ADataAcquisitionProtocolV2Checkpoint(unittest.TestCase):
    def test_happy_path(self) -> None:
        with tempfile.TemporaryDirectory() as td:
            td_path = Path(td)
            p1289 = td_path / "p1289.json"
            out = td_path / "p1290.json"
            p1289.write_text(
                json.dumps(
                    {
                        "next_priority": "R8A_DATA_ACQUISITION_PROTOCOL_V2",
                        "r8_program": {"recommended_first_move": "path_A_new_data"},
                    }
                ),
                encoding="utf-8",
            )
            subprocess.run(["python3", str(SCRIPT), "--p1289", str(p1289), "--out", str(out)], check=True)
            payload = json.loads(out.read_text(encoding="utf-8"))
            self.assertEqual(payload["r8a_protocol"]["status"], "DECLARED")
            self.assertEqual(payload["next_priority"], "R8A_EXECUTION_AND_MARGIN_REEVALUATION")

    def test_requires_path_a(self) -> None:
        with tempfile.TemporaryDirectory() as td:
            td_path = Path(td)
            p1289 = td_path / "p1289.json"
            out = td_path / "p1290.json"
            p1289.write_text(
                json.dumps(
                    {
                        "next_priority": "R8A_DATA_ACQUISITION_PROTOCOL_V2",
                        "r8_program": {"recommended_first_move": "path_B_formal_strengthening"},
                    }
                ),
                encoding="utf-8",
            )
            proc = subprocess.run(["python3", str(SCRIPT), "--p1289", str(p1289), "--out", str(out)], capture_output=True, text=True)
            self.assertNotEqual(proc.returncode, 0)
            self.assertIn("path_A_new_data", proc.stderr)


if __name__ == "__main__":
    unittest.main()
