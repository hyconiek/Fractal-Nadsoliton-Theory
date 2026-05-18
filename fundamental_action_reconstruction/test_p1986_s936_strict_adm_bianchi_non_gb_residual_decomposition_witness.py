from __future__ import annotations
import json, subprocess, unittest
from pathlib import Path

SCRIPT = Path(__file__).resolve().parent / "p1986_s936_strict_adm_bianchi_non_gb_residual_decomposition_witness.py"
OUT = Path(__file__).resolve().parent / "generated" / "p1986_s936_strict_adm_bianchi_non_gb_residual_decomposition_witness.json"

class TestP1986S936(unittest.TestCase):
    def test_decomposition_witness(self)->None:
        subprocess.run(["python3", str(SCRIPT)], check=True)
        p=json.loads(OUT.read_text(encoding="utf-8"))
        self.assertEqual(p["packet_id"], "P1986")
        self.assertEqual(p["result_kind"], "PASS_NON_GB_DECOMPOSITION_OBSTRUCTION_WITNESS")
        self.assertTrue(p["gatekeeper_checks"]["gb_channel_zero"])
        self.assertTrue(p["gatekeeper_checks"]["non_gb_channel_nonzero"])
        self.assertTrue(p["gatekeeper_checks"]["isotropic_limit_zero"])

if __name__=="__main__":
    unittest.main()
