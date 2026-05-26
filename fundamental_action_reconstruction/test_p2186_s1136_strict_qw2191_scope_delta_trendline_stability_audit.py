from __future__ import annotations

import json
import subprocess
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent
G = ROOT / "generated"


class TestP2186S1136StrictQW2191ScopeDeltaTrendlineStabilityAudit(unittest.TestCase):
    def test_packet_schema_and_trendline_export(self) -> None:
        subprocess.run(
            [sys.executable, str(ROOT / "p2186_s1136_strict_qw2191_scope_delta_trendline_stability_audit.py")],
            check=True,
        )
        data = json.loads(
            (
                G / "p2186_s1136_strict_qw2191_scope_delta_trendline_stability_audit.json"
            ).read_text(encoding="utf-8")
        )
        self.assertEqual(data["schema_version"], "p2186_s1136_v1")
        self.assertTrue(data["gatekeeper_checks"]["trendline_audit_exported"])
        self.assertIn("strict_qw2191_scope_delta_trendline_stability_audit", data)
        trendline = data["strict_qw2191_scope_delta_trendline_stability_audit"]["trendline"]
        self.assertIn("regression_count_series", trendline)


if __name__ == "__main__":
    unittest.main()
