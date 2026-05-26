from __future__ import annotations

import json
import subprocess
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent
G = ROOT / "generated"


class TestP2175S1125StrictQW2191StrongestEvidencePack(unittest.TestCase):
    def test_strongest_evidence_pack(self) -> None:
        subprocess.run(
            [sys.executable, str(ROOT / "p2175_s1125_strict_qw2191_missing_obligation_strongest_evidence_pack.py")],
            check=True,
        )
        d = json.loads(
            (G / "p2175_s1125_strict_qw2191_missing_obligation_strongest_evidence_pack.json").read_text(encoding="utf-8")
        )
        self.assertEqual(d["schema_version"], "p2175_s1125_v1")
        self.assertTrue(d["gatekeeper_checks"]["strongest_evidence_pack_exported"])
        pack = d["strict_qw2191_missing_obligation_strongest_evidence_pack"]
        self.assertIn("strongest_evidence_artifacts_for_passed", pack)
        self.assertIn("missing_obligation_remediation_map", pack)


if __name__ == "__main__":
    unittest.main()
