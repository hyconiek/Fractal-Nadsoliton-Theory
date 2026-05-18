from __future__ import annotations
import json, subprocess, unittest
from pathlib import Path

SCRIPT = Path(__file__).resolve().parent / "p2006_s956_strict_cutkosky_gx_backend_tensor_surrogate_and_covariance_classifier_witness.py"
OUT = Path(__file__).resolve().parent / "generated" / "p2006_s956_strict_cutkosky_gx_backend_tensor_surrogate_and_covariance_classifier_witness.json"

class TestP2006S956(unittest.TestCase):
    def test_tensor_covariance_classifier(self) -> None:
        subprocess.run(["python3", str(SCRIPT)], check=True)
        p = json.loads(OUT.read_text(encoding="utf-8"))
        self.assertEqual(p["packet_id"], "P2006")
        self.assertEqual(p["result_kind"], "PASS_GX_TENSOR_COVARIANCE_CLASSIFIER_WITNESS")
        self.assertTrue(p["gatekeeper_checks"]["tensor_surrogate_exported"])
        self.assertTrue(p["gatekeeper_checks"]["covariance_band_valid"])

if __name__ == "__main__":
    unittest.main()
