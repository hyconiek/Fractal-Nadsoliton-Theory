from __future__ import annotations
import json, subprocess, sys, unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent
G = ROOT / "generated"


class TestP2242S1192StrictNuBranchGroupPolicyEnforceableFloorConstraintProbe(unittest.TestCase):
    def test_packet(self) -> None:
        subprocess.run(
            [sys.executable, str(ROOT / "p2242_s1192_strict_nu_branch_group_policy_enforceable_floor_constraint_probe.py")],
            check=True,
        )
        data = json.loads(
            (G / "p2242_s1192_strict_nu_branch_group_policy_enforceable_floor_constraint_probe.json").read_text(
                encoding="utf-8"
            )
        )
        self.assertEqual(data["schema_version"], "p2242_s1192_v1")
        g = data["gatekeeper_checks"]
        self.assertTrue(g["enforceable_floor_constraint_exported"])
        self.assertTrue(g["mass_condition_computable"])
        self.assertTrue(g["uniform_floor_translation_exported"])


if __name__ == "__main__":
    unittest.main()
