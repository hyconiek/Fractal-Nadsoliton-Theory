from __future__ import annotations
import json, subprocess, unittest
from pathlib import Path
SCRIPT=Path(__file__).resolve().parent/"p2011_s961_strict_cutkosky_loop_amplitude_tensor_placeholder_and_coupled_scan_witness.py"
OUT=Path(__file__).resolve().parent/"generated"/"p2011_s961_strict_cutkosky_loop_amplitude_tensor_placeholder_and_coupled_scan_witness.json"
class TestP2011S961(unittest.TestCase):
    def test_placeholder(self)->None:
        subprocess.run(["python3",str(SCRIPT)],check=True)
        p=json.loads(OUT.read_text(encoding='utf-8'))
        self.assertEqual(p["packet_id"],"P2011")
        self.assertEqual(p["result_kind"],"PASS_LOOP_TENSOR_PLACEHOLDER_COUPLED_SCAN_WITNESS")
        self.assertTrue(p["gatekeeper_checks"]["loop_amplitude_tensor_placeholders_exported"])
if __name__=="__main__":
    unittest.main()
