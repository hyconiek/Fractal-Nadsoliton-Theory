from __future__ import annotations
import json, subprocess, unittest
from pathlib import Path

SCRIPT=Path(__file__).resolve().parent/"p1992_s942_strict_cutkosky_graviton_gauge_exact_phase_space_witness.py"
OUT=Path(__file__).resolve().parent/"generated"/"p1992_s942_strict_cutkosky_graviton_gauge_exact_phase_space_witness.json"

class TestP1992S942(unittest.TestCase):
    def test_phase_space_witness(self)->None:
        subprocess.run(["python3",str(SCRIPT)],check=True)
        p=json.loads(OUT.read_text(encoding="utf-8"))
        self.assertEqual(p["packet_id"],"P1992")
        self.assertEqual(p["result_kind"],"PASS_EXACT_PHASE_SPACE_BOUNDED_UNITARITY_PROXY_WITNESS")
        self.assertTrue(p["gatekeeper_checks"]["positive_imaginary_part_band_on_grid"])
        self.assertTrue(p["gatekeeper_checks"]["uncertainty_band_nonzero_and_bounded"])

if __name__=="__main__":
    unittest.main()
