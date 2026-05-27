from __future__ import annotations
import json, subprocess, sys, unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent
G = ROOT / "generated"


class TestP2270S1220StrictNuBranchGroupPolicySymbolicDerivativeBoundProbe(unittest.TestCase):
    def test_packet(self) -> None:
        subprocess.run(
            [sys.executable, str(ROOT / "p2270_s1220_strict_nu_branch_group_policy_symbolic_derivative_bound_probe.py")],
            check=True,
        )
        data = json.loads(
            (G / "p2270_s1220_strict_nu_branch_group_policy_symbolic_derivative_bound_probe.json").read_text(encoding="utf-8")
        )
        self.assertEqual(data["schema_version"], "p2270_s1220_v1")
        g = data["gatekeeper_checks"]
        self.assertTrue(g["symbolic_derivative_bounds_exported"])
        self.assertTrue(g["symbolic_monotone_in_rho_over_box"])
        self.assertTrue(g["symbolic_monotone_in_kappa_over_box"])
        self.assertTrue(g["symbolic_lipschitz_nonnegative"])


if __name__ == "__main__":
    unittest.main()
