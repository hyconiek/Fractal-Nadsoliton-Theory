from __future__ import annotations
import json, subprocess, unittest
from pathlib import Path

SCRIPT=Path(__file__).resolve().parent/"p1993_s943_strict_cutkosky_optical_defect_bounded_table_witness.py"
OUT=Path(__file__).resolve().parent/"generated"/"p1993_s943_strict_cutkosky_optical_defect_bounded_table_witness.json"

class TestP1993S943(unittest.TestCase):
    def test_bounded_defect(self)->None:
        subprocess.run(["python3",str(SCRIPT)],check=True)
        p=json.loads(OUT.read_text(encoding="utf-8"))
        self.assertEqual(p["packet_id"],"P1993")
        self.assertEqual(p["result_kind"],"PASS_OPTICAL_DEFECT_BOUNDED_TABLE_WITNESS")
        self.assertTrue(p["gatekeeper_checks"]["all_intervals_contain_zero"])
        self.assertTrue(p["gatekeeper_checks"]["proxy_residue_positive_on_selected_poles"])

if __name__=="__main__":
    unittest.main()
