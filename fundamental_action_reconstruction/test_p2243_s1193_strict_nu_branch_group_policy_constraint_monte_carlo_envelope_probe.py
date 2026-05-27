from __future__ import annotations
import json, subprocess, sys, unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent
G = ROOT / "generated"


class TestP2243S1193StrictNuBranchGroupPolicyConstraintMonteCarloEnvelopeProbe(unittest.TestCase):
    def test_packet(self) -> None:
        subprocess.run(
            [sys.executable, str(ROOT / "p2243_s1193_strict_nu_branch_group_policy_constraint_monte_carlo_envelope_probe.py")],
            check=True,
        )
        data = json.loads(
            (G / "p2243_s1193_strict_nu_branch_group_policy_constraint_monte_carlo_envelope_probe.json").read_text(
                encoding="utf-8"
            )
        )
        self.assertEqual(data["schema_version"], "p2243_s1193_v1")
        g = data["gatekeeper_checks"]
        self.assertTrue(g["monte_carlo_envelope_exported"])
        self.assertTrue(g["draw_count_sufficient"])
        self.assertTrue(g["violation_probability_computable"])


if __name__ == "__main__":
    unittest.main()
