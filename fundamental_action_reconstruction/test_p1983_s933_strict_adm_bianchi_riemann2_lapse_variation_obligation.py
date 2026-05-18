from __future__ import annotations

import json
import subprocess
import unittest
from pathlib import Path

SCRIPT = Path(__file__).resolve().parent / "p1983_s933_strict_adm_bianchi_riemann2_lapse_variation_obligation.py"
OUT = Path(__file__).resolve().parent / "generated" / "p1983_s933_strict_adm_bianchi_riemann2_lapse_variation_obligation.json"


class TestP1983S933(unittest.TestCase):
    def test_riemann2_lapse_variation_exports_higher_derivative_obligation(self) -> None:
        subprocess.run(["python3", str(SCRIPT)], check=True)
        payload = json.loads(OUT.read_text(encoding="utf-8"))

        self.assertEqual(payload["packet_id"], "P1983")
        self.assertEqual(payload["stage_id"], "S933")
        self.assertEqual(payload["result_kind"], "PASS_RIEMANN2_LAPSE_VARIATION_OBLIGATION_EXPORTED")
        self.assertEqual(payload["status"], "OPEN_OBSTRUCTION_WITH_TRACE")
        self.assertTrue(payload["depends_on"]["p1982_ricci2_lapse_obligation_present"])
        self.assertEqual(payload["riemann_block_setup"]["Kretschmann_formula"], "4*sum_i E_i^2 + 4*sum_{i<j} F_ij^2")
        self.assertTrue(payload["riemann2_lapse_euler_operator"]["isotropic_limit_zero"])
        self.assertTrue(payload["riemann2_lapse_euler_operator"]["contains_shear_acceleration"])
        self.assertTrue(payload["riemann2_lapse_euler_operator"]["contains_lapse_second_derivative"])
        self.assertTrue(payload["gatekeeper_checks"]["riemann2_lapse_obligation_passed"])
        self.assertIn("ToE closure", payload["theorem_export"]["not_licensed"])


if __name__ == "__main__":
    unittest.main()
