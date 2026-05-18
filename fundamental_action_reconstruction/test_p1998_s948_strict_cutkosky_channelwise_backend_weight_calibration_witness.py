from __future__ import annotations
import json, subprocess, unittest
from pathlib import Path

SCRIPT = Path(__file__).resolve().parent / "p1998_s948_strict_cutkosky_channelwise_backend_weight_calibration_witness.py"
OUT = Path(__file__).resolve().parent / "generated" / "p1998_s948_strict_cutkosky_channelwise_backend_weight_calibration_witness.json"


class TestP1998S948(unittest.TestCase):
    def test_backend_weight_calibration(self) -> None:
        subprocess.run(["python3", str(SCRIPT)], check=True)
        p = json.loads(OUT.read_text(encoding="utf-8"))
        self.assertEqual(p["packet_id"], "P1998")
        self.assertEqual(p["result_kind"], "PASS_BACKEND_WEIGHT_CALIBRATION_WITNESS")
        self.assertTrue(p["gatekeeper_checks"]["backend_weight_calibration"])
        self.assertTrue(p["gatekeeper_checks"]["weights_normalized"])


if __name__ == "__main__":
    unittest.main()
