from __future__ import annotations

import json
import subprocess
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent
G = ROOT / "generated"


class TestP2184S1134StrictQW2191ExternalAuditReplayDashboardScopeDeltaRegressionGuard(unittest.TestCase):
    def test_dashboard_scope_delta_regression_guard(self) -> None:
        subprocess.run(
            [
                sys.executable,
                str(ROOT / "p2184_s1134_strict_qw2191_external_audit_replay_dashboard_scope_delta_regression_guard.py"),
            ],
            check=True,
        )
        d = json.loads(
            (
                G
                / "p2184_s1134_strict_qw2191_external_audit_replay_dashboard_scope_delta_regression_guard.json"
            ).read_text(encoding="utf-8")
        )
        self.assertEqual(d["schema_version"], "p2184_s1134_v1")
        self.assertTrue(d["gatekeeper_checks"]["dashboard_exported"])
        self.assertIn(
            "strict_qw2191_external_audit_replay_dashboard_scope_delta_regression_guard",
            d,
        )


if __name__ == "__main__":
    unittest.main()
