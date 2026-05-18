from __future__ import annotations
import json, subprocess, unittest
from pathlib import Path
SCRIPT=Path(__file__).resolve().parent/"p2012_s962_strict_cutkosky_channel_loop_tensor_consistency_audit_witness.py"
OUT=Path(__file__).resolve().parent/"generated"/"p2012_s962_strict_cutkosky_channel_loop_tensor_consistency_audit_witness.json"
class TestP2012S962(unittest.TestCase):
    def test_audit(self)->None:
        subprocess.run(["python3",str(SCRIPT)],check=True)
        p=json.loads(OUT.read_text(encoding='utf-8'))
        self.assertEqual(p["packet_id"],"P2012")
        self.assertEqual(p["result_kind"],"PASS_CHANNEL_LOOP_TENSOR_CONSISTENCY_AUDIT_WITNESS")
        self.assertTrue(p["gatekeeper_checks"]["covariance_shapes_match"])
if __name__=="__main__":
    unittest.main()
