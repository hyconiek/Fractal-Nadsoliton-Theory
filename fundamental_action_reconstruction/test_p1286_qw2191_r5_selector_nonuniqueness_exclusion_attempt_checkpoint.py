from __future__ import annotations

import json
import subprocess
import tempfile
import unittest
from pathlib import Path

SCRIPT = Path(__file__).resolve().parent / "p1286_qw2191_r5_selector_nonuniqueness_exclusion_attempt_checkpoint.py"


class TestP1286QW2191R5SelectorNonuniquenessExclusionAttemptCheckpoint(unittest.TestCase):
    def test_happy_path(self) -> None:
        with tempfile.TemporaryDirectory() as td:
            td_path = Path(td)
            p1285 = td_path / "p1285.json"
            out = td_path / "p1286.json"
            p1285.write_text(
                json.dumps(
                    {
                        "next_priority": "R5_SELECTOR_NONUNIQUENESS_EXCLUSION_ATTEMPT",
                        "selector_source_identification": {"minimal_source_family": ["SSEL_SRC_A", "SSEL_SRC_B"]},
                    }
                ),
                encoding="utf-8",
            )
            subprocess.run(["python3", str(SCRIPT), "--p1285", str(p1285), "--out", str(out)], check=True)
            payload = json.loads(out.read_text(encoding="utf-8"))
            self.assertFalse(payload["r5_exclusion_attempt"]["nonuniqueness_excluded"])
            self.assertEqual(payload["next_priority"], "R6_NEW_OBSERVABLE_PROVIDER_OR_AXIOM_TAGGED_EXTENSION")

    def test_requires_two_classes(self) -> None:
        with tempfile.TemporaryDirectory() as td:
            td_path = Path(td)
            p1285 = td_path / "p1285.json"
            out = td_path / "p1286.json"
            p1285.write_text(
                json.dumps(
                    {
                        "next_priority": "R5_SELECTOR_NONUNIQUENESS_EXCLUSION_ATTEMPT",
                        "selector_source_identification": {"minimal_source_family": ["SSEL_SRC_A"]},
                    }
                ),
                encoding="utf-8",
            )
            proc = subprocess.run(["python3", str(SCRIPT), "--p1285", str(p1285), "--out", str(out)], capture_output=True, text=True)
            self.assertNotEqual(proc.returncode, 0)
            self.assertIn("at least two selector-source classes", proc.stderr)


if __name__ == "__main__":
    unittest.main()
