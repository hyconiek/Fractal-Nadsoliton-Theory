from __future__ import annotations

import json
import subprocess
import tempfile
import unittest
from pathlib import Path

SCRIPT = Path(__file__).resolve().parent / "p1304_qw2191_r20_nb1_nontransfer_lemma_proof_attempt_batch_1_checkpoint.py"


class TestP1304QW2191R20NB1NontransferLemmaProofAttemptBatch1Checkpoint(unittest.TestCase):
    def test_happy_path(self) -> None:
        with tempfile.TemporaryDirectory() as td:
            td_path = Path(td)
            p1303 = td_path / "p1303.json"
            out = td_path / "p1304.json"
            p1303.write_text(
                json.dumps(
                    {
                        "next_priority": "R20_NB1_NONTRANSFER_LEMMA_PROOF_ATTEMPT_BATCH_1",
                        "lemma_pack_status": "V1_DRAFT_COMPLETE",
                    }
                ),
                encoding="utf-8",
            )
            subprocess.run(["python3", str(SCRIPT), "--p1303", str(p1303), "--out", str(out)], check=True)
            payload = json.loads(out.read_text(encoding="utf-8"))
            self.assertEqual(payload["batch_status"], "PARTIAL_PROGRESS")

    def test_requires_v1_draft_complete(self) -> None:
        with tempfile.TemporaryDirectory() as td:
            td_path = Path(td)
            p1303 = td_path / "p1303.json"
            out = td_path / "p1304.json"
            p1303.write_text(
                json.dumps(
                    {
                        "next_priority": "R20_NB1_NONTRANSFER_LEMMA_PROOF_ATTEMPT_BATCH_1",
                        "lemma_pack_status": "BLOCKED",
                    }
                ),
                encoding="utf-8",
            )
            proc = subprocess.run(["python3", str(SCRIPT), "--p1303", str(p1303), "--out", str(out)], capture_output=True, text=True)
            self.assertNotEqual(proc.returncode, 0)
            self.assertIn("V1_DRAFT_COMPLETE", proc.stderr)


if __name__ == "__main__":
    unittest.main()
