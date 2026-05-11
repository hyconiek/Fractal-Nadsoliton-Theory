from __future__ import annotations

import json
import subprocess
import tempfile
import unittest
from pathlib import Path

SCRIPT = Path(__file__).resolve().parent / "p1305_qw2191_r21_nb1_lemma_gap_closure_provider_and_invariant_packet_checkpoint.py"


class TestP1305QW2191R21NB1LemmaGapClosureProviderAndInvariantPacketCheckpoint(unittest.TestCase):
    def test_happy_path(self) -> None:
        with tempfile.TemporaryDirectory() as td:
            td_path = Path(td)
            p1304 = td_path / "p1304.json"
            out = td_path / "p1305.json"
            p1304.write_text(
                json.dumps(
                    {
                        "next_priority": "R21_NB1_LEMMA_GAP_CLOSURE_PROVIDER_AND_INVARIANT_PACKET",
                        "batch_status": "PARTIAL_PROGRESS",
                    }
                ),
                encoding="utf-8",
            )
            subprocess.run(["python3", str(SCRIPT), "--p1304", str(p1304), "--out", str(out)], check=True)
            payload = json.loads(out.read_text(encoding="utf-8"))
            self.assertEqual(payload["r21_gap_closure_packet"]["status"], "PACKET_READY_FOR_BATCH2")

    def test_requires_partial_progress(self) -> None:
        with tempfile.TemporaryDirectory() as td:
            td_path = Path(td)
            p1304 = td_path / "p1304.json"
            out = td_path / "p1305.json"
            p1304.write_text(
                json.dumps(
                    {
                        "next_priority": "R21_NB1_LEMMA_GAP_CLOSURE_PROVIDER_AND_INVARIANT_PACKET",
                        "batch_status": "BLOCKED",
                    }
                ),
                encoding="utf-8",
            )
            proc = subprocess.run(["python3", str(SCRIPT), "--p1304", str(p1304), "--out", str(out)], capture_output=True, text=True)
            self.assertNotEqual(proc.returncode, 0)
            self.assertIn("PARTIAL_PROGRESS", proc.stderr)


if __name__ == "__main__":
    unittest.main()
