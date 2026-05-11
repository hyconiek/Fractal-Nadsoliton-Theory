from __future__ import annotations

import json
import subprocess
import tempfile
import unittest
from pathlib import Path

SCRIPT = Path(__file__).resolve().parent / "p1281_qw2191_r2_o1_common_gauge_transport_checkpoint.py"


class TestP1281QW2191R2O1CommonGaugeTransportCheckpoint(unittest.TestCase):
    def test_happy_path(self) -> None:
        with tempfile.TemporaryDirectory() as td:
            td_path = Path(td)
            p1280 = td_path / "p1280.json"
            out = td_path / "p1281.json"
            p1280.write_text(
                json.dumps(
                    {
                        "next_priority": "R2_O1_COMMON_GAUGE_TRANSPORT",
                        "closure_policy": {"global_qw2191_closure_allowed": False},
                    }
                ),
                encoding="utf-8",
            )
            subprocess.run(["python3", str(SCRIPT), "--p1280", str(p1280), "--out", str(out)], check=True)
            payload = json.loads(out.read_text(encoding="utf-8"))
            self.assertEqual(payload["packet"], "P1281")
            self.assertEqual(payload["r2_o1"]["status"], "PARTIAL_DISCHARGE")
            self.assertEqual(payload["next_priority"], "R2_O2_MISMATCH_CONTROL_LEMMA")
            self.assertFalse(payload["closure_policy"]["global_qw2191_closure_allowed"])

    def test_requires_priority_alignment(self) -> None:
        with tempfile.TemporaryDirectory() as td:
            td_path = Path(td)
            p1280 = td_path / "p1280.json"
            out = td_path / "p1281.json"
            p1280.write_text(
                json.dumps(
                    {
                        "next_priority": "OTHER",
                        "closure_policy": {"global_qw2191_closure_allowed": False},
                    }
                ),
                encoding="utf-8",
            )
            proc = subprocess.run(
                ["python3", str(SCRIPT), "--p1280", str(p1280), "--out", str(out)],
                capture_output=True,
                text=True,
            )
            self.assertNotEqual(proc.returncode, 0)
            self.assertIn("R2_O1_COMMON_GAUGE_TRANSPORT", proc.stderr)


if __name__ == "__main__":
    unittest.main()
