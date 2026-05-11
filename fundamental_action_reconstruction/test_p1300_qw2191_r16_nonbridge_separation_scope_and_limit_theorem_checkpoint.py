from __future__ import annotations

import json
import subprocess
import tempfile
import unittest
from pathlib import Path

SCRIPT = Path(__file__).resolve().parent / "p1300_qw2191_r16_nonbridge_separation_scope_and_limit_theorem_checkpoint.py"


class TestP1300QW2191R16NonbridgeSeparationScopeAndLimitTheoremCheckpoint(unittest.TestCase):
    def test_happy_path(self) -> None:
        with tempfile.TemporaryDirectory() as td:
            td_path = Path(td)
            p1299 = td_path / "p1299.json"
            out = td_path / "p1300.json"
            p1299.write_text(
                json.dumps(
                    {
                        "next_priority": "R16_NONBRIDGE_SEPARATION_SCOPE_AND_LIMIT_THEOREM",
                        "policy": {
                            "preferred_resolution_path": "NB1_NONBRIDGE",
                            "bridge_to_legacy_allowed": False,
                        },
                    }
                ),
                encoding="utf-8",
            )
            subprocess.run(["python3", str(SCRIPT), "--p1299", str(p1299), "--out", str(out)], check=True)
            payload = json.loads(out.read_text(encoding="utf-8"))
            self.assertEqual(payload["r16_nonbridge_scope_limit_theorem"]["status"], "SCOPE_LIMIT_DRAFTED")

    def test_requires_nb_path(self) -> None:
        with tempfile.TemporaryDirectory() as td:
            td_path = Path(td)
            p1299 = td_path / "p1299.json"
            out = td_path / "p1300.json"
            p1299.write_text(
                json.dumps(
                    {
                        "next_priority": "R16_NONBRIDGE_SEPARATION_SCOPE_AND_LIMIT_THEOREM",
                        "policy": {
                            "preferred_resolution_path": "B1_BRIDGE",
                            "bridge_to_legacy_allowed": False,
                        },
                    }
                ),
                encoding="utf-8",
            )
            proc = subprocess.run(["python3", str(SCRIPT), "--p1299", str(p1299), "--out", str(out)], capture_output=True, text=True)
            self.assertNotEqual(proc.returncode, 0)
            self.assertIn("NB1_NONBRIDGE", proc.stderr)


if __name__ == "__main__":
    unittest.main()
