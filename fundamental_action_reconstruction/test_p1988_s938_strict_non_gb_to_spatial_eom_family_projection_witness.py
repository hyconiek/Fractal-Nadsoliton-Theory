from __future__ import annotations
import json, subprocess, unittest
from pathlib import Path

SCRIPT = Path(__file__).resolve().parent / "p1988_s938_strict_non_gb_to_spatial_eom_family_projection_witness.py"
OUT = Path(__file__).resolve().parent / "generated" / "p1988_s938_strict_non_gb_to_spatial_eom_family_projection_witness.json"

class TestP1988S938(unittest.TestCase):
    def test_projection_witness(self) -> None:
        subprocess.run(["python3", str(SCRIPT)], check=True)
        p = json.loads(OUT.read_text(encoding="utf-8"))
        self.assertEqual(p["packet_id"], "P1988")
        self.assertEqual(p["result_kind"], "PASS_NON_GB_SPATIAL_PROJECTION_OBSTRUCTION_WITNESS")
        self.assertTrue(p["gatekeeper_checks"]["outside_eh_family_capacity_detected"])
        self.assertTrue(p["gatekeeper_checks"]["family_gap_norm_positive"])
        self.assertTrue(p["gatekeeper_checks"]["p1974_eh_residual_has_no_d2sigma_or_Ndd"])

if __name__ == "__main__":
    unittest.main()
