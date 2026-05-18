from __future__ import annotations
import json, subprocess, unittest
from pathlib import Path

SCRIPT = Path(__file__).resolve().parent / "p2007_s957_strict_cutkosky_tensor_vs_scalar_band_comparator_witness.py"
OUT = Path(__file__).resolve().parent / "generated" / "p2007_s957_strict_cutkosky_tensor_vs_scalar_band_comparator_witness.json"

class TestP2007S957(unittest.TestCase):
    def test_comparator(self) -> None:
        subprocess.run(["python3", str(SCRIPT)], check=True)
        p = json.loads(OUT.read_text(encoding="utf-8"))
        self.assertEqual(p["packet_id"], "P2007")
        self.assertEqual(p["result_kind"], "PASS_TENSOR_VS_SCALAR_COMPARATOR_WITNESS")
        self.assertTrue(p["gatekeeper_checks"]["scalar_scan_nonempty"])
        self.assertTrue(p["gatekeeper_checks"]["tensor_scan_nonempty"])

if __name__ == "__main__":
    unittest.main()
