from __future__ import annotations
import json, subprocess, sys, unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent
G = ROOT / "generated"


class TestP2246S1196StrictNuBranchGroupPolicyFeasibilitySurfaceProbe(unittest.TestCase):
    def test_packet(self) -> None:
        subprocess.run(
            [sys.executable, str(ROOT / "p2246_s1196_strict_nu_branch_group_policy_feasibility_surface_probe.py")],
            check=True,
        )
        data = json.loads(
            (G / "p2246_s1196_strict_nu_branch_group_policy_feasibility_surface_probe.json").read_text(encoding="utf-8")
        )
        self.assertEqual(data["schema_version"], "p2246_s1196_v1")
        g = data["gatekeeper_checks"]
        self.assertTrue(g["feasibility_surface_exported"])
        self.assertTrue(g["nonempty_surface_grid"])
        self.assertTrue(g["admissible_fraction_computable"])


if __name__ == "__main__":
    unittest.main()
