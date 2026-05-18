from __future__ import annotations
import json, subprocess, unittest
from pathlib import Path

SCRIPT=Path(__file__).resolve().parent/"p1997_s947_strict_cutkosky_channelwise_statesum_solver_witness.py"
OUT=Path(__file__).resolve().parent/"generated"/"p1997_s947_strict_cutkosky_channelwise_statesum_solver_witness.json"

class TestP1997S947(unittest.TestCase):
    def test_channelwise_statesum(self)->None:
        subprocess.run(["python3",str(SCRIPT)],check=True)
        p=json.loads(OUT.read_text(encoding="utf-8"))
        self.assertEqual(p["packet_id"],"P1997")
        self.assertEqual(p["result_kind"],"PASS_CHANNELWISE_STATESUM_DELTA_OPT_WITNESS")
        self.assertTrue(p["gatekeeper_checks"]["channelwise_statesum_exported"])
        self.assertTrue(p["gatekeeper_checks"]["nonzero_delta_detected"])

if __name__=="__main__":
    unittest.main()
