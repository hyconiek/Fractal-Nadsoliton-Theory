from __future__ import annotations

import json
import subprocess
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent
G = ROOT / "generated"


class TestP2177S1127StrictQW2191ObligationSpecificStrongestEvidenceAndMissingRecheck(unittest.TestCase):
    def test_obligation_specific_upgrade(self) -> None:
        subprocess.run(
            [sys.executable, str(ROOT / "p2177_s1127_strict_qw2191_obligation_specific_strongest_evidence_and_missing_recheck.py")],
            check=True,
        )
        d = json.loads(
            (G / "p2177_s1127_strict_qw2191_obligation_specific_strongest_evidence_and_missing_recheck.json").read_text(
                encoding="utf-8"
            )
        )
        self.assertEqual(d["schema_version"], "p2177_s1127_v1")
        self.assertTrue(d["gatekeeper_checks"]["obligation_specific_upgrade_exported"])
        mr = d["strict_qw2191_obligation_specific_strongest_evidence_and_missing_recheck"]["missing_recheck"]
        self.assertIn("remaining_missing_count", mr)


if __name__ == "__main__":
    unittest.main()
