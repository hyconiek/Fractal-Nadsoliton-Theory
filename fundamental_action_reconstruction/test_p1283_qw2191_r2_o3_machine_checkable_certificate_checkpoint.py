from __future__ import annotations

import json
import subprocess
import tempfile
import unittest
from pathlib import Path

SCRIPT = Path(__file__).resolve().parent / "p1283_qw2191_r2_o3_machine_checkable_certificate_checkpoint.py"


class TestP1283QW2191R2O3MachineCheckableCertificateCheckpoint(unittest.TestCase):
    def test_happy_path(self) -> None:
        with tempfile.TemporaryDirectory() as td:
            td_path = Path(td)
            p1282 = td_path / "p1282.json"
            out = td_path / "p1283.json"
            p1282.write_text(
                json.dumps(
                    {
                        "next_priority": "R2_O3_MACHINE_CHECKABLE_CERTIFICATE",
                        "evidence": {"epsilon_budget": {"eps_total_upper_bound": 0.025}},
                    }
                ),
                encoding="utf-8",
            )
            subprocess.run(["python3", str(SCRIPT), "--p1282", str(p1282), "--out", str(out)], check=True)
            payload = json.loads(out.read_text(encoding="utf-8"))
            self.assertEqual(payload["r2_status"], "O3_DISCHARGED")
            self.assertTrue(payload["certificate"]["all_checks_pass"])
            self.assertEqual(payload["next_priority"], "R3_INDEPENDENT_AUDIT_AND_REPLICATION")

    def test_requires_priority(self) -> None:
        with tempfile.TemporaryDirectory() as td:
            td_path = Path(td)
            p1282 = td_path / "p1282.json"
            out = td_path / "p1283.json"
            p1282.write_text(json.dumps({"next_priority": "OTHER", "evidence": {"epsilon_budget": {"eps_total_upper_bound": 0.025}}}), encoding="utf-8")
            proc = subprocess.run(["python3", str(SCRIPT), "--p1282", str(p1282), "--out", str(out)], capture_output=True, text=True)
            self.assertNotEqual(proc.returncode, 0)
            self.assertIn("R2_O3_MACHINE_CHECKABLE_CERTIFICATE", proc.stderr)


if __name__ == "__main__":
    unittest.main()
