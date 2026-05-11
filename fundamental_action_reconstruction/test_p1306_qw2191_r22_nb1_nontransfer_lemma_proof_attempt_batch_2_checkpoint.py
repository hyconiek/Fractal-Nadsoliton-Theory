from __future__ import annotations

import json
import subprocess
import tempfile
import unittest
from pathlib import Path

SCRIPT = Path(__file__).resolve().parent / "p1306_qw2191_r22_nb1_nontransfer_lemma_proof_attempt_batch_2_checkpoint.py"


class TestP1306QW2191R22NB1NontransferLemmaProofAttemptBatch2Checkpoint(unittest.TestCase):
    def test_happy_path(self) -> None:
        with tempfile.TemporaryDirectory() as td:
            td_path = Path(td)
            p1305 = td_path / "p1305.json"
            out = td_path / "p1306.json"
            p1305.write_text(
                json.dumps(
                    {
                        "next_priority": "R22_NB1_NONTRANSFER_LEMMA_PROOF_ATTEMPT_BATCH_2",
                        "r21_gap_closure_packet": {"status": "PACKET_READY_FOR_BATCH2"},
                    }
                ),
                encoding="utf-8",
            )
            subprocess.run(["python3", str(SCRIPT), "--p1305", str(p1305), "--out", str(out)], check=True)
            payload = json.loads(out.read_text(encoding="utf-8"))
            self.assertEqual(payload["batch_status"], "NEAR_COMPLETE")

    def test_requires_packet_ready(self) -> None:
        with tempfile.TemporaryDirectory() as td:
            td_path = Path(td)
            p1305 = td_path / "p1305.json"
            out = td_path / "p1306.json"
            p1305.write_text(
                json.dumps(
                    {
                        "next_priority": "R22_NB1_NONTRANSFER_LEMMA_PROOF_ATTEMPT_BATCH_2",
                        "r21_gap_closure_packet": {"status": "BLOCKED"},
                    }
                ),
                encoding="utf-8",
            )
            proc = subprocess.run(["python3", str(SCRIPT), "--p1305", str(p1305), "--out", str(out)], capture_output=True, text=True)
            self.assertNotEqual(proc.returncode, 0)
            self.assertIn("PACKET_READY_FOR_BATCH2", proc.stderr)


if __name__ == "__main__":
    unittest.main()
