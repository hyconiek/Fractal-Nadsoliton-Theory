from __future__ import annotations

import json
import subprocess
import tempfile
import unittest
from pathlib import Path

SCRIPT = Path(__file__).resolve().parent / "p1308_qw2191_r24_nb1_exotic_class_sweep_and_final_lemma_status_checkpoint.py"


class TestP1308QW2191R24NB1ExoticClassSweepAndFinalLemmaStatusCheckpoint(unittest.TestCase):
    def test_happy_path(self) -> None:
        with tempfile.TemporaryDirectory() as td:
            td_path = Path(td)
            p1307 = td_path / "p1307.json"
            out = td_path / "p1308.json"
            p1307.write_text(
                json.dumps(
                    {
                        "next_priority": "R24_NB1_EXOTIC_CLASS_SWEEP_AND_FINAL_LEMMA_STATUS",
                        "r23_conversion_protocol": {"status": "PROTOCOL_READY"},
                    }
                ),
                encoding="utf-8",
            )
            subprocess.run(["python3", str(SCRIPT), "--p1307", str(p1307), "--out", str(out)], check=True)
            payload = json.loads(out.read_text(encoding="utf-8"))
            self.assertEqual(payload["lnb1_2_final_status"], "PASS_STRICT")

    def test_requires_protocol_ready(self) -> None:
        with tempfile.TemporaryDirectory() as td:
            td_path = Path(td)
            p1307 = td_path / "p1307.json"
            out = td_path / "p1308.json"
            p1307.write_text(
                json.dumps(
                    {
                        "next_priority": "R24_NB1_EXOTIC_CLASS_SWEEP_AND_FINAL_LEMMA_STATUS",
                        "r23_conversion_protocol": {"status": "BLOCKED"},
                    }
                ),
                encoding="utf-8",
            )
            proc = subprocess.run(["python3", str(SCRIPT), "--p1307", str(p1307), "--out", str(out)], capture_output=True, text=True)
            self.assertNotEqual(proc.returncode, 0)
            self.assertIn("PROTOCOL_READY", proc.stderr)


if __name__ == "__main__":
    unittest.main()
