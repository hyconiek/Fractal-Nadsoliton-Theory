from __future__ import annotations

import json
import subprocess
import tempfile
import unittest
from pathlib import Path

SCRIPT = Path(__file__).resolve().parent / "p1297_qw2191_r13_governance_theorem_interface_preparation_checkpoint.py"


class TestP1297QW2191R13GovernanceTheoremInterfacePreparationCheckpoint(unittest.TestCase):
    def test_happy_path(self) -> None:
        with tempfile.TemporaryDirectory() as td:
            td_path = Path(td)
            p1296 = td_path / "p1296.json"
            out = td_path / "p1297.json"
            p1296.write_text(
                json.dumps(
                    {
                        "next_priority": "R13_GOVERNANCE_THEOREM_INTERFACE_PREPARATION",
                        "r12_motion_review": {"review_result": "CONDITIONAL_HOLD"},
                    }
                ),
                encoding="utf-8",
            )
            subprocess.run(["python3", str(SCRIPT), "--p1296", str(p1296), "--out", str(out)], check=True)
            payload = json.loads(out.read_text(encoding="utf-8"))
            self.assertEqual(payload["r13_interface"]["status"], "DRAFT_READY")
            self.assertEqual(payload["next_priority"], "R14_B1_NB1_INTERFACE_THEOREM_DRAFT")

    def test_requires_conditional_hold(self) -> None:
        with tempfile.TemporaryDirectory() as td:
            td_path = Path(td)
            p1296 = td_path / "p1296.json"
            out = td_path / "p1297.json"
            p1296.write_text(
                json.dumps(
                    {
                        "next_priority": "R13_GOVERNANCE_THEOREM_INTERFACE_PREPARATION",
                        "r12_motion_review": {"review_result": "APPROVED"},
                    }
                ),
                encoding="utf-8",
            )
            proc = subprocess.run(["python3", str(SCRIPT), "--p1296", str(p1296), "--out", str(out)], capture_output=True, text=True)
            self.assertNotEqual(proc.returncode, 0)
            self.assertIn("CONDITIONAL_HOLD", proc.stderr)


if __name__ == "__main__":
    unittest.main()
