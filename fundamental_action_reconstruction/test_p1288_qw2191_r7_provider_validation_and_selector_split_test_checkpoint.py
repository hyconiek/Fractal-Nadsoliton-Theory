from __future__ import annotations

import json
import subprocess
import tempfile
import unittest
from pathlib import Path

SCRIPT = Path(__file__).resolve().parent / "p1288_qw2191_r7_provider_validation_and_selector_split_test_checkpoint.py"


class TestP1288QW2191R7ProviderValidationAndSelectorSplitTestCheckpoint(unittest.TestCase):
    def test_happy_path(self) -> None:
        with tempfile.TemporaryDirectory() as td:
            td_path = Path(td)
            p1287 = td_path / "p1287.json"
            out = td_path / "p1288.json"
            p1287.write_text(
                json.dumps(
                    {
                        "next_priority": "R7_PROVIDER_VALIDATION_AND_SELECTOR_SPLIT_TEST",
                        "r6_program": {"observable_provider": {"status": "DECLARED", "id": "OBS_DELTA_PHASE_CURVATURE_V1"}},
                    }
                ),
                encoding="utf-8",
            )
            subprocess.run(["python3", str(SCRIPT), "--p1287", str(p1287), "--out", str(out)], check=True)
            payload = json.loads(out.read_text(encoding="utf-8"))
            self.assertEqual(payload["r7"]["provider_validation"]["calibration_result"], "PASS")
            self.assertEqual(payload["r7"]["selector_split_test"]["result"], "INCONCLUSIVE")

    def test_requires_declared_provider(self) -> None:
        with tempfile.TemporaryDirectory() as td:
            td_path = Path(td)
            p1287 = td_path / "p1287.json"
            out = td_path / "p1288.json"
            p1287.write_text(
                json.dumps({"next_priority": "R7_PROVIDER_VALIDATION_AND_SELECTOR_SPLIT_TEST", "r6_program": {"observable_provider": {"status": "MISSING"}}}),
                encoding="utf-8",
            )
            proc = subprocess.run(["python3", str(SCRIPT), "--p1287", str(p1287), "--out", str(out)], capture_output=True, text=True)
            self.assertNotEqual(proc.returncode, 0)
            self.assertIn("declared observable provider", proc.stderr)


if __name__ == "__main__":
    unittest.main()
