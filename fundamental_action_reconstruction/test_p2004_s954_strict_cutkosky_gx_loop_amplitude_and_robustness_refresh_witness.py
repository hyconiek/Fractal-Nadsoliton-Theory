from __future__ import annotations
import json, subprocess, unittest
from pathlib import Path

SCRIPT = Path(__file__).resolve().parent / "p2004_s954_strict_cutkosky_gx_loop_amplitude_and_robustness_refresh_witness.py"
OUT = Path(__file__).resolve().parent / "generated" / "p2004_s954_strict_cutkosky_gx_loop_amplitude_and_robustness_refresh_witness.json"


class TestP2004S954(unittest.TestCase):
    def test_gx_loop_refresh(self) -> None:
        subprocess.run(["python3", str(SCRIPT)], check=True)
        p = json.loads(OUT.read_text(encoding="utf-8"))
        self.assertEqual(p["packet_id"], "P2004")
        self.assertEqual(p["result_kind"], "PASS_GX_LOOP_AND_ROBUSTNESS_REFRESH_WITNESS")
        self.assertTrue(p["gatekeeper_checks"]["gx_loop_kernel_upgraded"])
        self.assertTrue(p["gatekeeper_checks"]["robustness_scan_nonempty"])


if __name__ == "__main__":
    unittest.main()
