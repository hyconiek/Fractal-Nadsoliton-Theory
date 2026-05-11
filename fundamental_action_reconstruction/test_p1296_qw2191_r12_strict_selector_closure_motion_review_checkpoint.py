from __future__ import annotations

import json
import subprocess
import tempfile
import unittest
from pathlib import Path

SCRIPT = Path(__file__).resolve().parent / "p1296_qw2191_r12_strict_selector_closure_motion_review_checkpoint.py"


class TestP1296QW2191R12StrictSelectorClosureMotionReviewCheckpoint(unittest.TestCase):
    def test_happy_path(self) -> None:
        with tempfile.TemporaryDirectory() as td:
            td_path = Path(td)
            p1295 = td_path / "p1295.json"
            out = td_path / "p1296.json"
            p1295.write_text(
                json.dumps(
                    {
                        "next_priority": "R12_STRICT_SELECTOR_CLOSURE_MOTION_REVIEW",
                        "r11": {"peer_replay": {"independent_run_status": "PASS"}},
                    }
                ),
                encoding="utf-8",
            )
            subprocess.run(["python3", str(SCRIPT), "--p1295", str(p1295), "--out", str(out)], check=True)
            payload = json.loads(out.read_text(encoding="utf-8"))
            self.assertEqual(payload["r12_motion_review"]["review_result"], "CONDITIONAL_HOLD")
            self.assertEqual(payload["next_priority"], "R13_GOVERNANCE_THEOREM_INTERFACE_PREPARATION")

    def test_requires_peer_replay_pass(self) -> None:
        with tempfile.TemporaryDirectory() as td:
            td_path = Path(td)
            p1295 = td_path / "p1295.json"
            out = td_path / "p1296.json"
            p1295.write_text(
                json.dumps(
                    {
                        "next_priority": "R12_STRICT_SELECTOR_CLOSURE_MOTION_REVIEW",
                        "r11": {"peer_replay": {"independent_run_status": "FAIL"}},
                    }
                ),
                encoding="utf-8",
            )
            proc = subprocess.run(["python3", str(SCRIPT), "--p1295", str(p1295), "--out", str(out)], capture_output=True, text=True)
            self.assertNotEqual(proc.returncode, 0)
            self.assertIn("PASS peer replay", proc.stderr)


if __name__ == "__main__":
    unittest.main()
