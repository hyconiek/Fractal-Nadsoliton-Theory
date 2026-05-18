from __future__ import annotations
import json, subprocess, unittest
from pathlib import Path

SCRIPT = Path(__file__).resolve().parent / "p2005_s955_strict_cutkosky_gx_backend_amplitude_bound_and_classifier_witness.py"
OUT = Path(__file__).resolve().parent / "generated" / "p2005_s955_strict_cutkosky_gx_backend_amplitude_bound_and_classifier_witness.json"


class TestP2005S955(unittest.TestCase):
    def test_backend_bound_classifier(self) -> None:
        subprocess.run(["python3", str(SCRIPT)], check=True)
        p = json.loads(OUT.read_text(encoding="utf-8"))
        self.assertEqual(p["packet_id"], "P2005")
        self.assertEqual(p["result_kind"], "PASS_GX_BACKEND_BOUND_CLASSIFIER_WITNESS")
        self.assertTrue(p["gatekeeper_checks"]["backend_bound_exported"])
        self.assertTrue(p["gatekeeper_checks"]["scan_nonempty"])


if __name__ == "__main__":
    unittest.main()
