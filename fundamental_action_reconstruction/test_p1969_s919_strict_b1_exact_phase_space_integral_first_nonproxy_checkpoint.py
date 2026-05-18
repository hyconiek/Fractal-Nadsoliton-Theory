from __future__ import annotations

import json
import subprocess
import unittest
from pathlib import Path

SCRIPT = Path(__file__).resolve().parent / "p1969_s919_strict_b1_exact_phase_space_integral_first_nonproxy_checkpoint.py"
OUT = Path(__file__).resolve().parent / "generated" / "p1969_s919_strict_b1_exact_phase_space_integral_first_nonproxy_checkpoint.json"


class TestP1969S919(unittest.TestCase):
    def test_exact_projected_integral_certificate(self) -> None:
        subprocess.run(["python3", str(SCRIPT)], check=True)
        payload = json.loads(OUT.read_text(encoding="utf-8"))

        self.assertEqual(payload["packet_id"], "P1969")
        self.assertEqual(payload["stage_id"], "S919")
        self.assertEqual(payload["status"], "OPEN_OBSTRUCTION_WITH_TRACE")
        self.assertTrue(payload["symbolic_result"]["closed_form_residual_is_zero"])
        self.assertEqual(payload["symbolic_result"]["closed_form_residual"], "0")
        self.assertTrue(payload["numeric_integral_table"]["all_positive_under_integration_error_only"])
        self.assertGreater(payload["numeric_integral_table"]["minimum_lower_bound"], 0.0)
        self.assertIn("not theorem-grade Cutkosky closure", payload["false_pass_guard"])


if __name__ == "__main__":
    unittest.main()
