from __future__ import annotations
import json, subprocess, unittest
from pathlib import Path

SCRIPT = Path(__file__).resolve().parent / "p2000_s950_strict_cutkosky_loop_kernel_two_channel_witness.py"
OUT = Path(__file__).resolve().parent / "generated" / "p2000_s950_strict_cutkosky_loop_kernel_two_channel_witness.json"


class TestP2000S950(unittest.TestCase):
    def test_loop_kernel_two_channel(self) -> None:
        subprocess.run(["python3", str(SCRIPT)], check=True)
        p = json.loads(OUT.read_text(encoding="utf-8"))
        self.assertEqual(p["packet_id"], "P2000")
        self.assertEqual(p["result_kind"], "PASS_LOOP_KERNEL_TWO_CHANNEL_WITNESS")
        self.assertTrue(p["gatekeeper_checks"]["loop_kernel_gg_exported"])
        self.assertTrue(p["gatekeeper_checks"]["loop_kernel_gh_exported"])
        self.assertTrue(p["gatekeeper_checks"]["hh_residual_proxy_marked"])


if __name__ == "__main__":
    unittest.main()
