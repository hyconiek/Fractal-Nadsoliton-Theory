from __future__ import annotations

import json
import subprocess
import unittest
from pathlib import Path

SCRIPT = Path(__file__).resolve().parent / "p1984_s934_strict_adm_bianchi_gauss_bonnet_lapse_cancellation_witness.py"
OUT = Path(__file__).resolve().parent / "generated" / "p1984_s934_strict_adm_bianchi_gauss_bonnet_lapse_cancellation_witness.json"


class TestP1984S934(unittest.TestCase):
    def test_gauss_bonnet_lapse_euler_cancellation(self) -> None:
        subprocess.run(["python3", str(SCRIPT)], check=True)
        payload = json.loads(OUT.read_text(encoding="utf-8"))

        self.assertEqual(payload["packet_id"], "P1984")
        self.assertEqual(payload["stage_id"], "S934")
        self.assertEqual(payload["result_kind"], "PASS_GB_LAPSE_CANCELLATION_WITNESS")
        self.assertEqual(payload["status"], "OPEN_OBSTRUCTION_WITH_TRACE")
        self.assertTrue(payload["depends_on"]["p1981_r2_lapse_obligation_present"])
        self.assertTrue(payload["depends_on"]["p1982_ricci2_lapse_obligation_present"])
        self.assertTrue(payload["depends_on"]["p1983_riemann2_lapse_obligation_present"])
        self.assertEqual(payload["gb_lapse_euler_operator"]["EL_N_GB_difference"], "0")
        self.assertTrue(payload["gb_lapse_euler_operator"]["higher_derivative_terms_cancelled_in_lapse_EL"])
        self.assertTrue(payload["gatekeeper_checks"]["gb_lapse_cancellation_witness_passed"])
        self.assertIn("ToE closure", payload["theorem_export"]["not_licensed"])


if __name__ == "__main__":
    unittest.main()
