from __future__ import annotations
import json, subprocess, sys, unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent
G = ROOT / "generated"


class TestP2265S1215StrictNuBranchTask3ClosureGapMatrixProbe(unittest.TestCase):
    def test_packet(self) -> None:
        subprocess.run(
            [sys.executable, str(ROOT / "p2265_s1215_strict_nu_branch_task3_closure_gap_matrix_probe.py")],
            check=True,
        )
        data = json.loads(
            (G / "p2265_s1215_strict_nu_branch_task3_closure_gap_matrix_probe.json").read_text(encoding="utf-8")
        )
        self.assertEqual(data["schema_version"], "p2265_s1215_v1")
        g = data["gatekeeper_checks"]
        self.assertTrue(g["task3_gap_matrix_exported"])
        self.assertTrue(g["closure_score_bounded"])
        self.assertTrue(g["contains_open_gaps"])


if __name__ == "__main__":
    unittest.main()
