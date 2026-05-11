from __future__ import annotations

import json
import subprocess
import tempfile
import unittest
from pathlib import Path

SCRIPT = Path(__file__).resolve().parent / "p1280_qw2191_r2_bound_transport_checkpoint.py"


class TestP1280QW2191R2BoundTransportCheckpoint(unittest.TestCase):
    def test_happy_path(self) -> None:
        with tempfile.TemporaryDirectory() as td:
            td_path = Path(td)
            p1279 = td_path / "p1279.json"
            out = td_path / "p1280.json"
            p1279.write_text(
                json.dumps(
                    {
                        "theorem": {"status": "PARTIAL_DISCHARGE"},
                        "next_priority": "R2_BOUND_TRANSPORT",
                    }
                ),
                encoding="utf-8",
            )
            subprocess.run(["python3", str(SCRIPT), "--p1279", str(p1279), "--out", str(out)], check=True)
            payload = json.loads(out.read_text(encoding="utf-8"))
            self.assertEqual(payload["packet"], "P1280")
            self.assertEqual(payload["r2_program"]["entry_gate"], "OPEN")
            self.assertFalse(payload["closure_policy"]["global_qw2191_closure_allowed"])
            self.assertEqual(payload["next_priority"], "R2_O1_COMMON_GAUGE_TRANSPORT")

    def test_requires_partial_discharge(self) -> None:
        with tempfile.TemporaryDirectory() as td:
            td_path = Path(td)
            p1279 = td_path / "p1279.json"
            out = td_path / "p1280.json"
            p1279.write_text(
                json.dumps(
                    {
                        "theorem": {"status": "UNSET"},
                        "next_priority": "R2_BOUND_TRANSPORT",
                    }
                ),
                encoding="utf-8",
            )
            proc = subprocess.run(
                ["python3", str(SCRIPT), "--p1279", str(p1279), "--out", str(out)],
                capture_output=True,
                text=True,
            )
            self.assertNotEqual(proc.returncode, 0)
            self.assertIn("PARTIAL_DISCHARGE", proc.stderr)


if __name__ == "__main__":
    unittest.main()
