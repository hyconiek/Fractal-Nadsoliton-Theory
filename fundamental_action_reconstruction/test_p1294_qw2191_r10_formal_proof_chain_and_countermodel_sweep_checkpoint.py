from __future__ import annotations

import json
import subprocess
import tempfile
import unittest
from pathlib import Path

SCRIPT = Path(__file__).resolve().parent / "p1294_qw2191_r10_formal_proof_chain_and_countermodel_sweep_checkpoint.py"


class TestP1294QW2191R10FormalProofChainAndCountermodelSweepCheckpoint(unittest.TestCase):
    def test_happy_path(self) -> None:
        with tempfile.TemporaryDirectory() as td:
            td_path = Path(td)
            p1293 = td_path / "p1293.json"
            out = td_path / "p1294.json"
            p1293.write_text(
                json.dumps(
                    {
                        "next_priority": "R10_FORMAL_PROOF_CHAIN_AND_COUNTERMODEL_SWEEP",
                        "r9_theorem_draft": {"status": "DRAFT"},
                    }
                ),
                encoding="utf-8",
            )
            subprocess.run(["python3", str(SCRIPT), "--p1293", str(p1293), "--out", str(out)], check=True)
            payload = json.loads(out.read_text(encoding="utf-8"))
            self.assertEqual(payload["r10"]["countermodel_sweep"]["countermodels_found"], 0)
            self.assertEqual(payload["next_priority"], "R11_FORMAL_PROOF_COMPLETION_AND_PEER_REPLAY")

    def test_requires_draft_status(self) -> None:
        with tempfile.TemporaryDirectory() as td:
            td_path = Path(td)
            p1293 = td_path / "p1293.json"
            out = td_path / "p1294.json"
            p1293.write_text(
                json.dumps(
                    {
                        "next_priority": "R10_FORMAL_PROOF_CHAIN_AND_COUNTERMODEL_SWEEP",
                        "r9_theorem_draft": {"status": "OPEN"},
                    }
                ),
                encoding="utf-8",
            )
            proc = subprocess.run(["python3", str(SCRIPT), "--p1293", str(p1293), "--out", str(out)], capture_output=True, text=True)
            self.assertNotEqual(proc.returncode, 0)
            self.assertIn("status DRAFT", proc.stderr)


if __name__ == "__main__":
    unittest.main()
