from __future__ import annotations

import json
import subprocess
import tempfile
import unittest
from pathlib import Path

SCRIPT = Path(__file__).resolve().parent / "p1282_qw2191_r2_o2_mismatch_control_lemma_checkpoint.py"


class TestP1282QW2191R2O2MismatchControlLemmaCheckpoint(unittest.TestCase):
    def test_happy_path(self) -> None:
        with tempfile.TemporaryDirectory() as td:
            td_path = Path(td)
            p1281 = td_path / "p1281.json"
            out = td_path / "p1282.json"
            p1281.write_text(
                json.dumps(
                    {
                        "next_priority": "R2_O2_MISMATCH_CONTROL_LEMMA",
                        "r2_o1": {
                            "transport_maps": [
                                {"map_id": "tau_1"},
                                {"map_id": "tau_2"},
                                {"map_id": "tau_3"},
                            ]
                        },
                    }
                ),
                encoding="utf-8",
            )
            subprocess.run(["python3", str(SCRIPT), "--p1281", str(p1281), "--out", str(out)], check=True)
            payload = json.loads(out.read_text(encoding="utf-8"))
            self.assertEqual(payload["packet"], "P1282")
            self.assertFalse(payload["closure_policy"]["legacy_strict_bridge_required"])
            self.assertEqual(payload["next_priority"], "R2_O3_MACHINE_CHECKABLE_CERTIFICATE")

    def test_requires_priority(self) -> None:
        with tempfile.TemporaryDirectory() as td:
            td_path = Path(td)
            p1281 = td_path / "p1281.json"
            out = td_path / "p1282.json"
            p1281.write_text(json.dumps({"next_priority": "OTHER", "r2_o1": {"transport_maps": [{"map_id": "tau_1"}, {"map_id": "tau_2"}]}}), encoding="utf-8")
            proc = subprocess.run(["python3", str(SCRIPT), "--p1281", str(p1281), "--out", str(out)], capture_output=True, text=True)
            self.assertNotEqual(proc.returncode, 0)
            self.assertIn("R2_O2_MISMATCH_CONTROL_LEMMA", proc.stderr)


if __name__ == "__main__":
    unittest.main()
