from __future__ import annotations

import json
import subprocess
import unittest
from pathlib import Path

SCRIPT = Path(__file__).resolve().parent / "p1980_s930_strict_adm_bianchi_eh_lapse_variation_witness.py"
OUT = Path(__file__).resolve().parent / "generated" / "p1980_s930_strict_adm_bianchi_eh_lapse_variation_witness.json"


class TestP1980S930(unittest.TestCase):
    def test_eh_lapse_variation_derives_minus_q_shear(self) -> None:
        subprocess.run(["python3", str(SCRIPT)], check=True)
        payload = json.loads(OUT.read_text(encoding="utf-8"))

        self.assertEqual(payload["packet_id"], "P1980")
        self.assertEqual(payload["stage_id"], "S930")
        self.assertEqual(payload["result_kind"], "PASS_EH_LAPSE_SHEAR_TERM_DERIVED")
        self.assertEqual(payload["status"], "OPEN_OBSTRUCTION_WITH_TRACE")
        self.assertTrue(payload["depends_on"]["p1979_current_export_audit_present"])
        self.assertEqual(payload["symbolic_variation"]["normalized_lapse_constraint_v"], "v1*v2 + v1*v3 + v2*v3")
        self.assertTrue(payload["symbolic_variation"]["variation_matches_pair_sum"])
        self.assertEqual(payload["shear_specialization"]["Gnn_BianchiI_EH_constraint"], "3*H**2 - sigma1**2 - sigma1*sigma2 - sigma2**2")
        self.assertEqual(payload["shear_specialization"]["Bianchi_minus_FRW"], "-sigma1**2 - sigma1*sigma2 - sigma2**2")
        self.assertTrue(payload["gatekeeper_checks"]["eh_lapse_witness_passed"])
        self.assertIn("ToE closure", payload["theorem_export"]["not_licensed"])


if __name__ == "__main__":
    unittest.main()
