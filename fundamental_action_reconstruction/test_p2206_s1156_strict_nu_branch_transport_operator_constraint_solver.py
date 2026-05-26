from __future__ import annotations

import json
import subprocess
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent
G = ROOT / "generated"


class TestP2206S1156StrictNuBranchTransportOperatorConstraintSolver(unittest.TestCase):
    def test_packet(self) -> None:
        subprocess.run(
            [sys.executable, str(ROOT / "p2206_s1156_strict_nu_branch_transport_operator_constraint_solver.py")],
            check=True,
        )
        data = json.loads(
            (G / "p2206_s1156_strict_nu_branch_transport_operator_constraint_solver.json").read_text(encoding="utf-8")
        )
        self.assertEqual(data["schema_version"], "p2206_s1156_v1")
        self.assertTrue(data["gatekeeper_checks"]["constraint_solver_exported"])
        verdict = data["strict_nu_branch_transport_operator_constraint_solver"]["solver_verdict"]
        self.assertIn("best_m_on_grid", verdict)


if __name__ == "__main__":
    unittest.main()
