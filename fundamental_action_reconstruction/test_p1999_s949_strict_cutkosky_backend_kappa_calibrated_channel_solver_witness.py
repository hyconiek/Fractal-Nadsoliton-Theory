from __future__ import annotations
import json, subprocess, unittest
from pathlib import Path

SCRIPT = Path(__file__).resolve().parent / "p1999_s949_strict_cutkosky_backend_kappa_calibrated_channel_solver_witness.py"
OUT = Path(__file__).resolve().parent / "generated" / "p1999_s949_strict_cutkosky_backend_kappa_calibrated_channel_solver_witness.json"


class TestP1999S949(unittest.TestCase):
    def test_backend_kappa_calibrated_solver(self) -> None:
        subprocess.run(["python3", str(SCRIPT)], check=True)
        p = json.loads(OUT.read_text(encoding="utf-8"))
        self.assertEqual(p["packet_id"], "P1999")
        self.assertEqual(p["result_kind"], "PASS_BACKEND_KAPPA_CALIBRATED_CHANNEL_SOLVER_WITNESS")
        self.assertTrue(p["gatekeeper_checks"]["backend_weight_calibration"])
        self.assertTrue(p["gatekeeper_checks"]["backend_kappa_calibration"])
        self.assertTrue(p["gatekeeper_checks"]["kappas_normalized"])


if __name__ == "__main__":
    unittest.main()
