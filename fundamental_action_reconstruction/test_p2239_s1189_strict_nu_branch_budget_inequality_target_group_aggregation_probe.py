from __future__ import annotations
import json, subprocess, sys, unittest
from pathlib import Path
ROOT = Path(__file__).resolve().parent
G = ROOT / "generated"

class TestP2239S1189StrictNuBranchBudgetInequalityTargetGroupAggregationProbe(unittest.TestCase):
    def test_packet(self) -> None:
        subprocess.run([sys.executable, str(ROOT / "p2239_s1189_strict_nu_branch_budget_inequality_target_group_aggregation_probe.py")], check=True)
        data = json.loads((G / "p2239_s1189_strict_nu_branch_budget_inequality_target_group_aggregation_probe.json").read_text(encoding="utf-8"))
        self.assertEqual(data["schema_version"], "p2239_s1189_v1")
        g = data["gatekeeper_checks"]
        self.assertTrue(g["target_group_aggregation_exported"])
        self.assertTrue(g["global_margin_nonnegative"])
        self.assertTrue(g["groupwise_compatibility_holds"])

if __name__ == "__main__":
    unittest.main()
