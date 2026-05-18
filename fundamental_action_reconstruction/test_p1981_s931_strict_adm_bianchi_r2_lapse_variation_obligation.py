from __future__ import annotations

import json
import subprocess
import unittest
from pathlib import Path

SCRIPT = Path(__file__).resolve().parent / "p1981_s931_strict_adm_bianchi_r2_lapse_variation_obligation.py"
OUT = Path(__file__).resolve().parent / "generated" / "p1981_s931_strict_adm_bianchi_r2_lapse_variation_obligation.json"


class TestP1981S931(unittest.TestCase):
    def test_r2_lapse_variation_exports_nontrivial_obligation(self) -> None:
        subprocess.run(["python3", str(SCRIPT)], check=True)
        payload = json.loads(OUT.read_text(encoding="utf-8"))

        self.assertEqual(payload["packet_id"], "P1981")
        self.assertEqual(payload["stage_id"], "S931")
        self.assertEqual(payload["result_kind"], "PASS_R2_LAPSE_VARIATION_OBLIGATION_EXPORTED")
        self.assertEqual(payload["status"], "OPEN_OBSTRUCTION_WITH_TRACE")
        self.assertTrue(payload["depends_on"]["p1980_eh_lapse_witness_present"])
        self.assertEqual(payload["ricci_scalar_setup"]["R_BianchiI_minus_R_FRW"], "2*Q/N**2")
        self.assertEqual(payload["r2_lapse_euler_operator"]["EL_N_difference"], "-12*V*(6*H**2*Q - 2*H*Qd + 4*Hd*Q + Q**2)/N**4")
        self.assertTrue(payload["r2_lapse_euler_operator"]["contains_Qdot"])
        self.assertTrue(payload["r2_lapse_euler_operator"]["contains_Q_squared"])
        self.assertTrue(payload["gatekeeper_checks"]["r2_lapse_obligation_passed"])
        self.assertIn("ToE closure", payload["theorem_export"]["not_licensed"])


if __name__ == "__main__":
    unittest.main()
