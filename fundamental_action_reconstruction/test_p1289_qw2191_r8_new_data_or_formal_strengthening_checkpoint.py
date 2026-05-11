from __future__ import annotations

import json
import subprocess
import tempfile
import unittest
from pathlib import Path

SCRIPT = Path(__file__).resolve().parent / "p1289_qw2191_r8_new_data_or_formal_strengthening_checkpoint.py"


class TestP1289QW2191R8NewDataOrFormalStrengtheningCheckpoint(unittest.TestCase):
    def test_happy_path(self) -> None:
        with tempfile.TemporaryDirectory() as td:
            td_path = Path(td)
            p1288 = td_path / "p1288.json"
            out = td_path / "p1289.json"
            p1288.write_text(
                json.dumps(
                    {
                        "next_priority": "R8_NEW_DATA_OR_FORMAL_STRENGTHENING",
                        "r7": {"selector_split_test": {"result": "INCONCLUSIVE"}},
                    }
                ),
                encoding="utf-8",
            )
            subprocess.run(["python3", str(SCRIPT), "--p1288", str(p1288), "--out", str(out)], check=True)
            payload = json.loads(out.read_text(encoding="utf-8"))
            self.assertEqual(payload["r8_program"]["recommended_first_move"], "path_A_new_data")
            self.assertEqual(payload["next_priority"], "R8A_DATA_ACQUISITION_PROTOCOL_V2")

    def test_requires_inconclusive_split(self) -> None:
        with tempfile.TemporaryDirectory() as td:
            td_path = Path(td)
            p1288 = td_path / "p1288.json"
            out = td_path / "p1289.json"
            p1288.write_text(
                json.dumps(
                    {
                        "next_priority": "R8_NEW_DATA_OR_FORMAL_STRENGTHENING",
                        "r7": {"selector_split_test": {"result": "PASS"}},
                    }
                ),
                encoding="utf-8",
            )
            proc = subprocess.run(["python3", str(SCRIPT), "--p1288", str(p1288), "--out", str(out)], capture_output=True, text=True)
            self.assertNotEqual(proc.returncode, 0)
            self.assertIn("INCONCLUSIVE", proc.stderr)


if __name__ == "__main__":
    unittest.main()
