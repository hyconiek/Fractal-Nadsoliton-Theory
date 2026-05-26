from __future__ import annotations

import json
import subprocess
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent
G = ROOT / "generated"


class TestP2173S1123StrictQW2191ObligationIndependentCrosschecker(unittest.TestCase):
    def test_crosschecker(self) -> None:
        subprocess.run([sys.executable, str(ROOT / "p2173_s1123_strict_qw2191_obligation_independent_crosschecker.py")], check=True)
        d = json.loads((G / "p2173_s1123_strict_qw2191_obligation_independent_crosschecker.json").read_text(encoding="utf-8"))
        self.assertEqual(d["schema_version"], "p2173_s1123_v1")
        self.assertTrue(d["gatekeeper_checks"]["crosscheck_exported"])
        self.assertTrue(d["gatekeeper_checks"]["all_consistent"])


if __name__ == "__main__":
    unittest.main()
