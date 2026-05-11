from __future__ import annotations

import json
import subprocess
import tempfile
import unittest
from pathlib import Path

SCRIPT = Path(__file__).resolve().parent / "p1293_qw2191_r9_formal_selector_source_theorem_draft_checkpoint.py"


class TestP1293QW2191R9FormalSelectorSourceTheoremDraftCheckpoint(unittest.TestCase):
    def test_happy_path(self) -> None:
        with tempfile.TemporaryDirectory() as td:
            td_path = Path(td)
            p1292 = td_path / "p1292.json"
            out = td_path / "p1293.json"
            p1292.write_text(
                json.dumps(
                    {
                        "next_priority": "R9_FORMAL_SELECTOR_SOURCE_THEOREM_DRAFT",
                        "r8b_split_rerun": {"report": {"result": "DECISIVE_FOR_SSEL_SRC_A"}},
                    }
                ),
                encoding="utf-8",
            )
            subprocess.run(["python3", str(SCRIPT), "--p1292", str(p1292), "--out", str(out)], check=True)
            payload = json.loads(out.read_text(encoding="utf-8"))
            self.assertEqual(payload["r9_theorem_draft"]["status"], "DRAFT")
            self.assertEqual(payload["next_priority"], "R10_FORMAL_PROOF_CHAIN_AND_COUNTERMODEL_SWEEP")

    def test_requires_decisive_result(self) -> None:
        with tempfile.TemporaryDirectory() as td:
            td_path = Path(td)
            p1292 = td_path / "p1292.json"
            out = td_path / "p1293.json"
            p1292.write_text(
                json.dumps(
                    {
                        "next_priority": "R9_FORMAL_SELECTOR_SOURCE_THEOREM_DRAFT",
                        "r8b_split_rerun": {"report": {"result": "INCONCLUSIVE"}},
                    }
                ),
                encoding="utf-8",
            )
            proc = subprocess.run(["python3", str(SCRIPT), "--p1292", str(p1292), "--out", str(out)], capture_output=True, text=True)
            self.assertNotEqual(proc.returncode, 0)
            self.assertIn("DECISIVE_FOR_SSEL_SRC_A", proc.stderr)


if __name__ == "__main__":
    unittest.main()
