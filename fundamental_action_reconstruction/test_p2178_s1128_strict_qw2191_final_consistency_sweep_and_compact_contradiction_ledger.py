from __future__ import annotations

import json
import subprocess
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent
G = ROOT / "generated"


class TestP2178S1128StrictQW2191FinalConsistencySweep(unittest.TestCase):
    def test_final_sweep(self) -> None:
        subprocess.run(
            [sys.executable, str(ROOT / "p2178_s1128_strict_qw2191_final_consistency_sweep_and_compact_contradiction_ledger.py")],
            check=True,
        )
        d = json.loads(
            (G / "p2178_s1128_strict_qw2191_final_consistency_sweep_and_compact_contradiction_ledger.json").read_text(
                encoding="utf-8"
            )
        )
        self.assertEqual(d["schema_version"], "p2178_s1128_v1")
        self.assertTrue(d["gatekeeper_checks"]["final_sweep_exported"])
        self.assertIn(
            "compact_ledger",
            d["strict_qw2191_final_consistency_sweep_and_compact_contradiction_ledger"],
        )


if __name__ == "__main__":
    unittest.main()
