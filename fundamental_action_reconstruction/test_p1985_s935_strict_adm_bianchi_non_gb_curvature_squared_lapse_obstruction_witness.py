from __future__ import annotations

import json
import subprocess
import unittest
from pathlib import Path

SCRIPT = Path(__file__).resolve().parent / "p1985_s935_strict_adm_bianchi_non_gb_curvature_squared_lapse_obstruction_witness.py"
OUT = Path(__file__).resolve().parent / "generated" / "p1985_s935_strict_adm_bianchi_non_gb_curvature_squared_lapse_obstruction_witness.json"


class TestP1985S935(unittest.TestCase):
    def test_non_gb_lapse_obstruction_witness(self) -> None:
        subprocess.run(["python3", str(SCRIPT)], check=True)
        payload = json.loads(OUT.read_text(encoding="utf-8"))

        self.assertEqual(payload["packet_id"], "P1985")
        self.assertEqual(payload["stage_id"], "S935")
        self.assertEqual(payload["result_kind"], "PASS_NON_GB_LAPSE_OBSTRUCTION_WITNESS")
        self.assertEqual(payload["status"], "OPEN_OBSTRUCTION_WITH_TRACE")
        self.assertTrue(payload["gatekeeper_checks"]["obstruction_witness_passed"])
        self.assertTrue(payload["weighted_non_gb_lapse_operator"]["contains_high_derivatives"])
        self.assertFalse(payload["weighted_non_gb_lapse_operator"]["is_identically_zero"])
        self.assertTrue(payload["weighted_non_gb_lapse_operator"]["isotropic_limit_zero"])
        self.assertIn("Q", payload["strict_symbolic_identifications"])
        self.assertIn("Qd", payload["strict_symbolic_identifications"])


if __name__ == "__main__":
    unittest.main()
