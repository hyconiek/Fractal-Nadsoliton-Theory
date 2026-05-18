from __future__ import annotations
import json, subprocess, unittest
from pathlib import Path

SCRIPT = Path(__file__).resolve().parent / "p2009_s959_strict_cutkosky_channelwise_tensor_coupled_covariance_classifier_witness.py"
OUT = Path(__file__).resolve().parent / "generated" / "p2009_s959_strict_cutkosky_channelwise_tensor_coupled_covariance_classifier_witness.json"

class TestP2009S959(unittest.TestCase):
    def test_coupled_covariance_classifier(self) -> None:
        subprocess.run(["python3", str(SCRIPT)], check=True)
        p = json.loads(OUT.read_text(encoding="utf-8"))
        self.assertEqual(p["packet_id"], "P2009")
        self.assertEqual(p["result_kind"], "PASS_CHANNELWISE_COUPLED_COVARIANCE_CLASSIFIER_WITNESS")
        self.assertTrue(p["gatekeeper_checks"]["channel_covariance_exported"])
        self.assertTrue(p["gatekeeper_checks"]["covariance_psd"])

if __name__ == "__main__":
    unittest.main()
