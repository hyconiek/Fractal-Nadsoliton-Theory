from __future__ import annotations
import json, subprocess, unittest
from pathlib import Path

SCRIPT=Path(__file__).resolve().parent/"p1995_s945_strict_cutkosky_scheme_sensitivity_and_optical_defect_witness.py"
OUT=Path(__file__).resolve().parent/"generated"/"p1995_s945_strict_cutkosky_scheme_sensitivity_and_optical_defect_witness.json"

class TestP1995S945(unittest.TestCase):
    def test_scheme_sensitivity(self)->None:
        subprocess.run(["python3",str(SCRIPT)],check=True)
        p=json.loads(OUT.read_text(encoding="utf-8"))
        self.assertEqual(p["packet_id"],"P1995")
        self.assertEqual(p["result_kind"],"PASS_SCHEME_SENSITIVITY_OPTICAL_DEFECT_WITNESS")
        self.assertTrue(p["gatekeeper_checks"]["nonzero_optical_defect_detected_under_scheme_perturbation"])
        self.assertTrue(p["gatekeeper_checks"]["residue_positive_on_scanned_family"])

if __name__=="__main__":
    unittest.main()
