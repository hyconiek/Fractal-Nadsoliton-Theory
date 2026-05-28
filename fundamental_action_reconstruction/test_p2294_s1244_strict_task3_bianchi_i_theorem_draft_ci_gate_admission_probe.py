from __future__ import annotations
import json
import subprocess
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent
G = ROOT / "generated"


class TestP2294S1244StrictTask3BianchiITheoremDraftCiGateAdmissionProbe(unittest.TestCase):
    def test_packet(self) -> None:
        subprocess.run(
            [sys.executable, str(ROOT / "p2294_s1244_strict_task3_bianchi_i_theorem_draft_ci_gate_admission_probe.py")],
            check=True,
        )
        data = json.loads(
            (G / "p2294_s1244_strict_task3_bianchi_i_theorem_draft_ci_gate_admission_probe.json").read_text(encoding="utf-8")
        )
        self.assertEqual(data["schema_version"], "p2294_s1244_v1")
        g = data["gatekeeper_checks"]
        self.assertTrue(g["ci_gate_decision_exported"])
        self.assertTrue(g["theorem_attempt_decision_exported"])
        self.assertTrue(g["open_implies_metadata_accept"])
        self.assertTrue(g["open_implies_negative_controls_ok"])
        probe = data["strict_task3_bianchi_i_theorem_draft_ci_gate_admission_probe"]
        self.assertEqual(probe["ci_gate_decision"], "CI_GATE_BLOCK")
        self.assertEqual(probe["theorem_attempt_decision"], "THEOREM_DRAFT_HOLD")


if __name__ == "__main__":
    unittest.main()
