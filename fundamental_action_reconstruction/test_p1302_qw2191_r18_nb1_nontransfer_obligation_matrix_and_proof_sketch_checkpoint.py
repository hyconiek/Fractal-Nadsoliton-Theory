from __future__ import annotations

import json
import subprocess
import tempfile
import unittest
from pathlib import Path

SCRIPT = Path(__file__).resolve().parent / "p1302_qw2191_r18_nb1_nontransfer_obligation_matrix_and_proof_sketch_checkpoint.py"


class TestP1302QW2191R18NB1NontransferObligationMatrixAndProofSketchCheckpoint(unittest.TestCase):
    def test_happy_path(self) -> None:
        with tempfile.TemporaryDirectory() as td:
            td_path = Path(td)
            p1301 = td_path / "p1301.json"
            out = td_path / "p1302.json"
            p1301.write_text(
                json.dumps(
                    {
                        "next_priority": "R18_NB1_NONTRANSFER_OBLIGATION_MATRIX_AND_PROOF_SKETCH",
                        "r17_nb1_nontransfer_theorem": {"status": "DRAFT_WITH_OBLIGATIONS"},
                    }
                ),
                encoding="utf-8",
            )
            subprocess.run(["python3", str(SCRIPT), "--p1301", str(p1301), "--out", str(out)], check=True)
            payload = json.loads(out.read_text(encoding="utf-8"))
            self.assertEqual(payload["proof_sketch_status"], "MATRIX_DRAFTED")
            self.assertEqual(len(payload["r18_obligation_matrix"]), 3)

    def test_requires_draft_with_obligations(self) -> None:
        with tempfile.TemporaryDirectory() as td:
            td_path = Path(td)
            p1301 = td_path / "p1301.json"
            out = td_path / "p1302.json"
            p1301.write_text(
                json.dumps(
                    {
                        "next_priority": "R18_NB1_NONTRANSFER_OBLIGATION_MATRIX_AND_PROOF_SKETCH",
                        "r17_nb1_nontransfer_theorem": {"status": "BLOCKED"},
                    }
                ),
                encoding="utf-8",
            )
            proc = subprocess.run(["python3", str(SCRIPT), "--p1301", str(p1301), "--out", str(out)], capture_output=True, text=True)
            self.assertNotEqual(proc.returncode, 0)
            self.assertIn("DRAFT_WITH_OBLIGATIONS", proc.stderr)


if __name__ == "__main__":
    unittest.main()
