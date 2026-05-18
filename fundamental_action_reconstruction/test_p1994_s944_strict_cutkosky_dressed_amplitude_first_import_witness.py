from __future__ import annotations
import json, subprocess, unittest
from pathlib import Path

SCRIPT=Path(__file__).resolve().parent/"p1994_s944_strict_cutkosky_dressed_amplitude_first_import_witness.py"
OUT=Path(__file__).resolve().parent/"generated"/"p1994_s944_strict_cutkosky_dressed_amplitude_first_import_witness.json"

class TestP1994S944(unittest.TestCase):
    def test_dressed_import(self)->None:
        subprocess.run(["python3",str(SCRIPT)],check=True)
        p=json.loads(OUT.read_text(encoding="utf-8"))
        self.assertEqual(p["packet_id"],"P1994")
        self.assertEqual(p["result_kind"],"PASS_DRESSED_IMPORT_OPTICAL_ZERO_PROXY_WITNESS")
        self.assertTrue(p["gatekeeper_checks"]["optical_defect_zero_on_grid"])
        self.assertTrue(p["gatekeeper_checks"]["proxy_residue_positive_on_poles"])

if __name__=="__main__":
    unittest.main()
