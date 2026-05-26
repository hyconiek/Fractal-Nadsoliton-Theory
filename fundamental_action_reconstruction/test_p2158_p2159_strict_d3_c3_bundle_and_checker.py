from __future__ import annotations

import json
import subprocess
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent
G = ROOT / "generated"


class TestP2158P2159StrictD3C3BundleAndChecker(unittest.TestCase):
    def test_export(self) -> None:
        subprocess.run([sys.executable, str(ROOT / "p2158_s1108_strict_d3_c3_transport_theorem_grade_derivation_bundle.py")], check=True)
        subprocess.run([sys.executable, str(ROOT / "p2159_s1109_strict_d3_c3_transport_independent_symbolic_checker.py")], check=True)

        d8 = json.loads((G / "p2158_s1108_strict_d3_c3_transport_theorem_grade_derivation_bundle.json").read_text(encoding="utf-8"))
        d9 = json.loads((G / "p2159_s1109_strict_d3_c3_transport_independent_symbolic_checker.json").read_text(encoding="utf-8"))

        self.assertEqual(d8["schema_version"], "p2158_s1108_v1")
        self.assertEqual(d9["schema_version"], "p2159_s1109_v1")
        self.assertTrue(d8["gatekeeper_checks"]["derivation_bundle_exported"])
        self.assertTrue(d9["gatekeeper_checks"]["independent_checker_exported"])
        self.assertTrue(d8["gatekeeper_checks"]["full_d3_covariance_transport_proven"])
        self.assertTrue(d8["gatekeeper_checks"]["c3_theorem_proven"])
        self.assertTrue(d9["gatekeeper_checks"]["full_d3_covariance_transport_proven"])
        self.assertTrue(d9["gatekeeper_checks"]["c3_theorem_proven"])


if __name__ == "__main__":
    unittest.main()
