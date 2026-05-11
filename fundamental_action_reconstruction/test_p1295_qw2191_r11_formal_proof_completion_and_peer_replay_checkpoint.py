from __future__ import annotations

import json
import subprocess
import tempfile
import unittest
from pathlib import Path

SCRIPT = Path(__file__).resolve().parent / "p1295_qw2191_r11_formal_proof_completion_and_peer_replay_checkpoint.py"


class TestP1295QW2191R11FormalProofCompletionAndPeerReplayCheckpoint(unittest.TestCase):
    def test_happy_path(self) -> None:
        with tempfile.TemporaryDirectory() as td:
            td_path = Path(td)
            p1294 = td_path / "p1294.json"
            out = td_path / "p1295.json"
            p1294.write_text(
                json.dumps(
                    {
                        "next_priority": "R11_FORMAL_PROOF_COMPLETION_AND_PEER_REPLAY",
                        "r10": {"countermodel_sweep": {"countermodels_found": 0}},
                    }
                ),
                encoding="utf-8",
            )
            subprocess.run(["python3", str(SCRIPT), "--p1294", str(p1294), "--out", str(out)], check=True)
            payload = json.loads(out.read_text(encoding="utf-8"))
            self.assertEqual(payload["r11"]["peer_replay"]["independent_run_status"], "PASS")
            self.assertEqual(payload["next_priority"], "R12_STRICT_SELECTOR_CLOSURE_MOTION_REVIEW")

    def test_requires_zero_countermodels(self) -> None:
        with tempfile.TemporaryDirectory() as td:
            td_path = Path(td)
            p1294 = td_path / "p1294.json"
            out = td_path / "p1295.json"
            p1294.write_text(
                json.dumps(
                    {
                        "next_priority": "R11_FORMAL_PROOF_COMPLETION_AND_PEER_REPLAY",
                        "r10": {"countermodel_sweep": {"countermodels_found": 1}},
                    }
                ),
                encoding="utf-8",
            )
            proc = subprocess.run(["python3", str(SCRIPT), "--p1294", str(p1294), "--out", str(out)], capture_output=True, text=True)
            self.assertNotEqual(proc.returncode, 0)
            self.assertIn("zero countermodels", proc.stderr)


if __name__ == "__main__":
    unittest.main()
