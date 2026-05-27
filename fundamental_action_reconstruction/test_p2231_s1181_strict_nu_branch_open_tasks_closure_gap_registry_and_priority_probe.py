from __future__ import annotations

import json
import subprocess
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent
G = ROOT / "generated"


class TestP2231S1181StrictNuBranchOpenTasksClosureGapRegistryAndPriorityProbe(unittest.TestCase):
    def test_packet(self) -> None:
        subprocess.run(
            [sys.executable, str(ROOT / "p2231_s1181_strict_nu_branch_open_tasks_closure_gap_registry_and_priority_probe.py")],
            check=True,
        )
        data = json.loads(
            (G / "p2231_s1181_strict_nu_branch_open_tasks_closure_gap_registry_and_priority_probe.json").read_text(encoding="utf-8")
        )
        self.assertEqual(data["schema_version"], "p2231_s1181_v1")
        g = data["gatekeeper_checks"]
        self.assertTrue(g["open_tasks_registry_exported"])
        self.assertTrue(g["priority_ranking_nonempty"])
        self.assertTrue(g["contains_global_transport_task"])


if __name__ == "__main__":
    unittest.main()
