from __future__ import annotations

import json
import subprocess
import tempfile
import unittest
from pathlib import Path

SCRIPT = Path(__file__).resolve().parent / "p1299_qw2191_r15_axiom_tagged_closure_policy_checkpoint.py"


class TestP1299QW2191R15AxiomTaggedClosurePolicyCheckpoint(unittest.TestCase):
    def test_happy_path(self) -> None:
        with tempfile.TemporaryDirectory() as td:
            td_path = Path(td)
            p1298 = td_path / "p1298.json"
            out = td_path / "p1299.json"
            p1298.write_text(
                json.dumps(
                    {
                        "next_priority": "R15_B1_NB1_OBLIGATION_MATRIX_AND_PROOF_PLAN",
                        "r14_interface_theorem": {"status": "THEOREM_INTERFACE_DRAFTED"},
                    }
                ),
                encoding="utf-8",
            )
            subprocess.run(["python3", str(SCRIPT), "--p1298", str(p1298), "--out", str(out)], check=True)
            payload = json.loads(out.read_text(encoding="utf-8"))
            self.assertFalse(payload["policy"]["bridge_to_legacy_allowed"])
            self.assertEqual(payload["policy"]["preferred_resolution_path"], "NB1_NONBRIDGE")
            self.assertEqual(payload["policy"]["axiom_augmented_closure_label"], "NON_STRICT")

    def test_requires_theorem_interface_drafted(self) -> None:
        with tempfile.TemporaryDirectory() as td:
            td_path = Path(td)
            p1298 = td_path / "p1298.json"
            out = td_path / "p1299.json"
            p1298.write_text(
                json.dumps(
                    {
                        "next_priority": "R15_B1_NB1_OBLIGATION_MATRIX_AND_PROOF_PLAN",
                        "r14_interface_theorem": {"status": "BLOCKED"},
                    }
                ),
                encoding="utf-8",
            )
            proc = subprocess.run(["python3", str(SCRIPT), "--p1298", str(p1298), "--out", str(out)], capture_output=True, text=True)
            self.assertNotEqual(proc.returncode, 0)
            self.assertIn("THEOREM_INTERFACE_DRAFTED", proc.stderr)


if __name__ == "__main__":
    unittest.main()
