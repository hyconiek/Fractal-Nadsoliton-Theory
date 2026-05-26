from __future__ import annotations

import json
import subprocess
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent
G = ROOT / "generated"


class TestP2185S1135StrictQW2191PeriodicCiSentinelAndNonbridgeScopeTagger(unittest.TestCase):
    def test_packet_schema_and_scope_tags(self) -> None:
        subprocess.run(
            [
                sys.executable,
                str(ROOT / "p2185_s1135_strict_qw2191_periodic_ci_sentinel_and_nonbridge_scope_tagger.py"),
            ],
            check=True,
        )
        data = json.loads(
            (
                G
                / "p2185_s1135_strict_qw2191_periodic_ci_sentinel_and_nonbridge_scope_tagger.json"
            ).read_text(encoding="utf-8")
        )
        self.assertEqual(data["schema_version"], "p2185_s1135_v1")
        self.assertTrue(data["gatekeeper_checks"]["ci_sentinel_exported"])
        self.assertTrue(data["gatekeeper_checks"]["nonbridge_scope_tags_applied"])
        self.assertIn("strict_qw2191_periodic_ci_sentinel_and_nonbridge_scope_tagger", data)


if __name__ == "__main__":
    unittest.main()
