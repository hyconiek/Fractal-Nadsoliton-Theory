from __future__ import annotations
import json, subprocess, unittest
from pathlib import Path

SCRIPT = Path(__file__).resolve().parent / "p2001_s951_strict_cutkosky_full_three_loop_kernel_plus_extra_channel_witness.py"
OUT = Path(__file__).resolve().parent / "generated" / "p2001_s951_strict_cutkosky_full_three_loop_kernel_plus_extra_channel_witness.json"


class TestP2001S951(unittest.TestCase):
    def test_full_three_plus_extra(self) -> None:
        subprocess.run(["python3", str(SCRIPT)], check=True)
        p = json.loads(OUT.read_text(encoding="utf-8"))
        self.assertEqual(p["packet_id"], "P2001")
        self.assertEqual(p["result_kind"], "PASS_FULL_THREE_LOOP_PLUS_EXTRA_CHANNEL_WITNESS")
        self.assertTrue(p["gatekeeper_checks"]["loop_kernel_hh_exported"])
        self.assertTrue(p["gatekeeper_checks"]["extra_channel_gx_exported"])


if __name__ == "__main__":
    unittest.main()
