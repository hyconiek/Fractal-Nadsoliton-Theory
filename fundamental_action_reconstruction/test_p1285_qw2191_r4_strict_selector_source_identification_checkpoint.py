from __future__ import annotations

import json
import subprocess
import tempfile
import unittest
from pathlib import Path

SCRIPT = Path(__file__).resolve().parent / "p1285_qw2191_r4_strict_selector_source_identification_checkpoint.py"


class TestP1285QW2191R4StrictSelectorSourceIdentificationCheckpoint(unittest.TestCase):
    def test_happy_path(self) -> None:
        with tempfile.TemporaryDirectory() as td:
            td_path = Path(td)
            p1284 = td_path / "p1284.json"
            out = td_path / "p1285.json"
            p1284.write_text(
                json.dumps(
                    {
                        "next_priority": "R4_STRICT_SELECTOR_SOURCE_IDENTIFICATION",
                        "independent_audit": {"result": "PASS"},
                    }
                ),
                encoding="utf-8",
            )
            subprocess.run(["python3", str(SCRIPT), "--p1284", str(p1284), "--out", str(out)], check=True)
            payload = json.loads(out.read_text(encoding="utf-8"))
            self.assertEqual(payload["selector_source_identification"]["status"], "PARTIAL_DISCHARGE")
            self.assertEqual(payload["next_priority"], "R5_SELECTOR_NONUNIQUENESS_EXCLUSION_ATTEMPT")

    def test_requires_audit_pass(self) -> None:
        with tempfile.TemporaryDirectory() as td:
            td_path = Path(td)
            p1284 = td_path / "p1284.json"
            out = td_path / "p1285.json"
            p1284.write_text(
                json.dumps(
                    {
                        "next_priority": "R4_STRICT_SELECTOR_SOURCE_IDENTIFICATION",
                        "independent_audit": {"result": "FAIL"},
                    }
                ),
                encoding="utf-8",
            )
            proc = subprocess.run(["python3", str(SCRIPT), "--p1284", str(p1284), "--out", str(out)], capture_output=True, text=True)
            self.assertNotEqual(proc.returncode, 0)
            self.assertIn("PASS independent audit", proc.stderr)


if __name__ == "__main__":
    unittest.main()
