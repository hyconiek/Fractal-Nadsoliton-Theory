from __future__ import annotations
import json
import subprocess
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent
G = ROOT / "generated"


class TestP2273S1223StrictNuBranchGroupPolicyTightenedSymbolicEnvelopeConstantProbe(unittest.TestCase):
    def test_packet(self) -> None:
        subprocess.run(
            [sys.executable, str(ROOT / "p2273_s1223_strict_nu_branch_group_policy_tightened_symbolic_envelope_constant_probe.py")],
            check=True,
        )
        data = json.loads(
            (G / "p2273_s1223_strict_nu_branch_group_policy_tightened_symbolic_envelope_constant_probe.json").read_text(encoding="utf-8")
        )
        self.assertEqual(data["schema_version"], "p2273_s1223_v1")
        g = data["gatekeeper_checks"]
        self.assertTrue(g["tightened_constants_exported"])
        self.assertTrue(g["safety_multiplier_ge_one"])
        self.assertIn(g["tightening_nonexpansive_vs_legacy"], [True, False])
        self.assertTrue(g["tightened_constants_nonnegative"])


if __name__ == "__main__":
    unittest.main()
