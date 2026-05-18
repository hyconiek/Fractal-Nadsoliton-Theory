from __future__ import annotations
import json, subprocess, unittest
from pathlib import Path

SCRIPT = Path(__file__).resolve().parent / "p1987_s937_strict_non_gb_residual_term_classification_witness.py"
OUT = Path(__file__).resolve().parent / "generated" / "p1987_s937_strict_non_gb_residual_term_classification_witness.json"

class TestP1987S937(unittest.TestCase):
    def test_term_classification(self) -> None:
        subprocess.run(["python3", str(SCRIPT)], check=True)
        p = json.loads(OUT.read_text(encoding="utf-8"))
        self.assertEqual(p["packet_id"], "P1987")
        self.assertEqual(p["result_kind"], "PASS_NON_GB_TERM_CLASSIFICATION_WITNESS")
        self.assertTrue(p["gatekeeper_checks"]["key_term_families_nonzero"])
        self.assertTrue(p["gatekeeper_checks"]["anisotropy_orders_2_and_4_present"])
        self.assertTrue(p["gatekeeper_checks"]["numeric_probe_nonzero"])

if __name__ == "__main__":
    unittest.main()
