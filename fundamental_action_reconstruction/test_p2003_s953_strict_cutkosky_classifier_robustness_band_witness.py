from __future__ import annotations
import json, subprocess, unittest
from pathlib import Path

SCRIPT = Path(__file__).resolve().parent / "p2003_s953_strict_cutkosky_classifier_robustness_band_witness.py"
OUT = Path(__file__).resolve().parent / "generated" / "p2003_s953_strict_cutkosky_classifier_robustness_band_witness.json"


class TestP2003S953(unittest.TestCase):
    def test_robustness_band(self) -> None:
        subprocess.run(["python3", str(SCRIPT)], check=True)
        p = json.loads(OUT.read_text(encoding="utf-8"))
        self.assertEqual(p["packet_id"], "P2003")
        self.assertEqual(p["result_kind"], "PASS_CLASSIFIER_ROBUSTNESS_BAND_WITNESS")
        self.assertTrue(p["gatekeeper_checks"]["nonempty_scan"])
        self.assertTrue(p["gatekeeper_checks"]["finite_band"])


if __name__ == "__main__":
    unittest.main()
