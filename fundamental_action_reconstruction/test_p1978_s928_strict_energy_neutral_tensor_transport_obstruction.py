from __future__ import annotations

import json
import subprocess
import unittest
from pathlib import Path

SCRIPT = Path(__file__).resolve().parent / "p1978_s928_strict_energy_neutral_tensor_transport_obstruction.py"
OUT = Path(__file__).resolve().parent / "generated" / "p1978_s928_strict_energy_neutral_tensor_transport_obstruction.json"


class TestP1978S928(unittest.TestCase):
    def test_energy_neutral_transport_obstruction(self) -> None:
        subprocess.run(["python3", str(SCRIPT)], check=True)
        payload = json.loads(OUT.read_text(encoding="utf-8"))

        self.assertEqual(payload["packet_id"], "P1978")
        self.assertEqual(payload["stage_id"], "S928")
        self.assertEqual(payload["result_kind"], "PASS_BOUNDED_OBSTRUCTION")
        self.assertEqual(payload["status"], "OPEN_OBSTRUCTION_WITH_TRACE")
        self.assertTrue(payload["depends_on"]["p1974_present"])
        self.assertTrue(payload["depends_on"]["p1977_bounded_no_go_present"])
        self.assertEqual(payload["symbolic_core"]["S_traceless_sum"], "0")
        self.assertEqual(payload["symbolic_core"]["energy_neutral_leftover_G00"], "-sigma1**2 - sigma1*sigma2 - sigma2**2")
        self.assertIn(payload["symbolic_core"]["tracefree_spatial_leftover_sum"], {"-3*sigma1**2 - 3*sigma1*sigma2 - 3*sigma2**2", "-3*(sigma1**2 + sigma1*sigma2 + sigma2**2)"})
        self.assertTrue(payload["gatekeeper_checks"]["energy_neutral_leaves_G00_minus_Q"])
        self.assertTrue(payload["gatekeeper_checks"]["tracefree_spatial_leaves_sum_minus_3Q"])
        self.assertTrue(payload["gatekeeper_checks"]["bounded_obstruction_passed"])
        self.assertIn("ToE closure", payload["theorem_export"]["not_licensed"])


if __name__ == "__main__":
    unittest.main()
