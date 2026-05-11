from __future__ import annotations

import json
import subprocess
import tempfile
import unittest
from pathlib import Path

SCRIPT = Path(__file__).resolve().parent / "p1303_qw2191_r19_nb1_nontransfer_lemma_pack_v1_checkpoint.py"


class TestP1303QW2191R19NB1NontransferLemmaPackV1Checkpoint(unittest.TestCase):
    def test_happy_path(self) -> None:
        with tempfile.TemporaryDirectory() as td:
            td_path = Path(td)
            p1302 = td_path / "p1302.json"
            out = td_path / "p1303.json"
            p1302.write_text(
                json.dumps(
                    {
                        "next_priority": "R19_NB1_NONTRANSFER_LEMMA_PACK_V1",
                        "proof_sketch_status": "MATRIX_DRAFTED",
                    }
                ),
                encoding="utf-8",
            )
            subprocess.run(["python3", str(SCRIPT), "--p1302", str(p1302), "--out", str(out)], check=True)
            payload = json.loads(out.read_text(encoding="utf-8"))
            self.assertEqual(payload["lemma_pack_status"], "V1_DRAFT_COMPLETE")

    def test_requires_matrix_drafted(self) -> None:
        with tempfile.TemporaryDirectory() as td:
            td_path = Path(td)
            p1302 = td_path / "p1302.json"
            out = td_path / "p1303.json"
            p1302.write_text(
                json.dumps(
                    {
                        "next_priority": "R19_NB1_NONTRANSFER_LEMMA_PACK_V1",
                        "proof_sketch_status": "BLOCKED",
                    }
                ),
                encoding="utf-8",
            )
            proc = subprocess.run(["python3", str(SCRIPT), "--p1302", str(p1302), "--out", str(out)], capture_output=True, text=True)
            self.assertNotEqual(proc.returncode, 0)
            self.assertIn("MATRIX_DRAFTED", proc.stderr)


if __name__ == "__main__":
    unittest.main()
