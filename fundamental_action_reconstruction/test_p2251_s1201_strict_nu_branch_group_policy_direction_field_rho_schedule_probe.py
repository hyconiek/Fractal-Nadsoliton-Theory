from __future__ import annotations
import json, subprocess, sys, unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent
G = ROOT / "generated"


class TestP2251S1201StrictNuBranchGroupPolicyDirectionFieldRhoScheduleProbe(unittest.TestCase):
    def test_packet(self) -> None:
        subprocess.run(
            [sys.executable, str(ROOT / "p2251_s1201_strict_nu_branch_group_policy_direction_field_rho_schedule_probe.py")],
            check=True,
        )
        data = json.loads(
            (G / "p2251_s1201_strict_nu_branch_group_policy_direction_field_rho_schedule_probe.json").read_text(
                encoding="utf-8"
            )
        )
        self.assertEqual(data["schema_version"], "p2251_s1201_v1")
        g = data["gatekeeper_checks"]
        self.assertTrue(g["scenario_comparison_exported"])
        self.assertTrue(g["all_scenarios_nonnegative_margin"])
        self.assertTrue(g["best_scenario_identified"])


if __name__ == "__main__":
    unittest.main()
