from __future__ import annotations
import json
import subprocess
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent
G = ROOT / "generated"


class TestP2293S1243StrictTask3BianchiIMetadataValidatorNegativeControlMatrixProbe(unittest.TestCase):
    def test_packet(self) -> None:
        subprocess.run(
            [sys.executable, str(ROOT / "p2293_s1243_strict_task3_bianchi_i_metadata_validator_negative_control_matrix_probe.py")],
            check=True,
        )
        data = json.loads(
            (G / "p2293_s1243_strict_task3_bianchi_i_metadata_validator_negative_control_matrix_probe.json").read_text(encoding="utf-8")
        )
        self.assertEqual(data["schema_version"], "p2293_s1243_v1")
        g = data["gatekeeper_checks"]
        self.assertTrue(g["matrix_rows_exported"])
        self.assertTrue(g["all_expectations_met"])
        self.assertTrue(g["contains_accept_case"])
        self.assertTrue(g["contains_reject_cases"])


if __name__ == "__main__":
    unittest.main()
