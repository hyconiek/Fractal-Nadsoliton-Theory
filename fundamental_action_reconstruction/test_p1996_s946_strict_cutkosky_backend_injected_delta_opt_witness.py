from __future__ import annotations
import json, subprocess, unittest
from pathlib import Path

SCRIPT=Path(__file__).resolve().parent/"p1996_s946_strict_cutkosky_backend_injected_delta_opt_witness.py"
OUT=Path(__file__).resolve().parent/"generated"/"p1996_s946_strict_cutkosky_backend_injected_delta_opt_witness.json"

class TestP1996S946(unittest.TestCase):
    def test_backend_injected(self)->None:
        subprocess.run(["python3",str(SCRIPT)],check=True)
        p=json.loads(OUT.read_text(encoding="utf-8"))
        self.assertEqual(p["packet_id"],"P1996")
        self.assertEqual(p["result_kind"],"PASS_BACKEND_INJECTED_DELTA_OPT_WITNESS")
        self.assertTrue(p["gatekeeper_checks"]["backend_profile_injected"])
        self.assertTrue(p["gatekeeper_checks"]["nonzero_delta_detected"])

if __name__=="__main__":
    unittest.main()
