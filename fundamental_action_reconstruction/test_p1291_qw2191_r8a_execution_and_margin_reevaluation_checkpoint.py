from __future__ import annotations

import json
import subprocess
import tempfile
import unittest
from pathlib import Path

SCRIPT = Path(__file__).resolve().parent / "p1291_qw2191_r8a_execution_and_margin_reevaluation_checkpoint.py"


class TestP1291QW2191R8AExecutionAndMarginReevaluationCheckpoint(unittest.TestCase):
    def test_happy_path(self) -> None:
        with tempfile.TemporaryDirectory() as td:
            td_path = Path(td)
            p1290 = td_path / "p1290.json"
            out = td_path / "p1291.json"
            p1290.write_text(
                json.dumps(
                    {
                        "next_priority": "R8A_EXECUTION_AND_MARGIN_REEVALUATION",
                        "r8a_protocol": {
                            "sample_size_min": 48,
                            "acceptance_criteria": {
                                "noise_floor_max": 0.006,
                                "margin_gain_min": 0.008,
                                "cross_sector_consistency_min": 0.95,
                            },
                        },
                    }
                ),
                encoding="utf-8",
            )
            subprocess.run(["python3", str(SCRIPT), "--p1290", str(p1290), "--out", str(out)], check=True)
            payload = json.loads(out.read_text(encoding="utf-8"))
            self.assertTrue(payload["r8a_execution"]["acceptance_pass"])
            self.assertEqual(payload["next_priority"], "R8B_SELECTOR_SPLIT_RERUN_WITH_UPDATED_MARGIN")

    def test_requires_priority(self) -> None:
        with tempfile.TemporaryDirectory() as td:
            td_path = Path(td)
            p1290 = td_path / "p1290.json"
            out = td_path / "p1291.json"
            p1290.write_text(json.dumps({"next_priority": "OTHER", "r8a_protocol": {"sample_size_min": 48}}), encoding="utf-8")
            proc = subprocess.run(["python3", str(SCRIPT), "--p1290", str(p1290), "--out", str(out)], capture_output=True, text=True)
            self.assertNotEqual(proc.returncode, 0)
            self.assertIn("R8A_EXECUTION_AND_MARGIN_REEVALUATION", proc.stderr)


if __name__ == "__main__":
    unittest.main()
