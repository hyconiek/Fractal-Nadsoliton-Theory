from __future__ import annotations
import json, subprocess, unittest
from pathlib import Path

SCRIPT = Path(__file__).resolve().parent / "p2008_s958_strict_cutkosky_direct_backend_tensor_object_classifier_witness.py"
OUT = Path(__file__).resolve().parent / "generated" / "p2008_s958_strict_cutkosky_direct_backend_tensor_object_classifier_witness.json"

class TestP2008S958(unittest.TestCase):
    def test_direct_tensor_classifier(self) -> None:
        subprocess.run(["python3", str(SCRIPT)], check=True)
        p = json.loads(OUT.read_text(encoding="utf-8"))
        self.assertEqual(p["packet_id"], "P2008")
        self.assertEqual(p["result_kind"], "PASS_DIRECT_BACKEND_TENSOR_CLASSIFIER_WITNESS")
        self.assertTrue(p["gatekeeper_checks"]["direct_backend_tensor_exported"])
        self.assertTrue(p["gatekeeper_checks"]["tensor_psd"])

if __name__ == "__main__":
    unittest.main()
