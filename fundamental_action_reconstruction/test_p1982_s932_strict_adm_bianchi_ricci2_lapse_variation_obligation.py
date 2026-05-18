from __future__ import annotations

import json
import subprocess
import unittest
from pathlib import Path

SCRIPT = Path(__file__).resolve().parent / "p1982_s932_strict_adm_bianchi_ricci2_lapse_variation_obligation.py"
OUT = Path(__file__).resolve().parent / "generated" / "p1982_s932_strict_adm_bianchi_ricci2_lapse_variation_obligation.json"


class TestP1982S932(unittest.TestCase):
    def test_ricci2_lapse_variation_exports_higher_derivative_obligation(self) -> None:
        subprocess.run(["python3", str(SCRIPT)], check=True)
        payload = json.loads(OUT.read_text(encoding="utf-8"))

        self.assertEqual(payload["packet_id"], "P1982")
        self.assertEqual(payload["stage_id"], "S932")
        self.assertEqual(payload["result_kind"], "PASS_RICCI2_LAPSE_VARIATION_OBLIGATION_EXPORTED")
        self.assertEqual(payload["status"], "OPEN_OBSTRUCTION_WITH_TRACE")
        self.assertTrue(payload["depends_on"]["p1981_r2_lapse_obligation_present"])
        self.assertTrue(payload["ricci2_lapse_euler_operator"]["isotropic_limit_zero"])
        self.assertTrue(payload["ricci2_lapse_euler_operator"]["contains_shear_acceleration"])
        self.assertTrue(payload["ricci2_lapse_euler_operator"]["contains_lapse_second_derivative"])
        self.assertTrue(payload["ricci2_lapse_euler_operator"]["contains_shear_velocity_squared"])
        self.assertTrue(payload["gatekeeper_checks"]["ricci2_lapse_obligation_passed"])
        self.assertIn("ToE closure", payload["theorem_export"]["not_licensed"])


if __name__ == "__main__":
    unittest.main()
