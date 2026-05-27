from __future__ import annotations
import json
import subprocess
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent
G = ROOT / "generated"


class TestP2290S1240StrictTask3BianchiIDeterministicTheoremPrecheckEvaluatorProbe(unittest.TestCase):
    def test_packet(self) -> None:
        subprocess.run(
            [sys.executable, str(ROOT / "p2290_s1240_strict_task3_bianchi_i_deterministic_theorem_precheck_evaluator_probe.py")],
            check=True,
        )
        data = json.loads(
            (G / "p2290_s1240_strict_task3_bianchi_i_deterministic_theorem_precheck_evaluator_probe.json").read_text(encoding="utf-8")
        )
        self.assertEqual(data["schema_version"], "p2290_s1240_v1")
        g = data["gatekeeper_checks"]
        self.assertTrue(g["checks_exported"])
        self.assertTrue(g["fingerprint_length_ok"])
        self.assertTrue(g["decision_valid"])


if __name__ == "__main__":
    unittest.main()
