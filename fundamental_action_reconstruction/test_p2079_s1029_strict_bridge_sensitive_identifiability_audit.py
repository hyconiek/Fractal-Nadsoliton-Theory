from __future__ import annotations

import json
import subprocess
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent
SCRIPTS = [
    ROOT / "p2078_s1028_strict_constrained_multiobjective_robustness_audit.py",
    ROOT / "p2079_s1029_strict_bridge_sensitive_identifiability_audit.py",
]
OUT = ROOT / "generated" / "p2079_s1029_strict_bridge_sensitive_identifiability_audit.json"


class TestP2079BridgeSensitiveIdentifiabilityAudit(unittest.TestCase):
    def test_p2079_exports_identifiability_and_keeps_c3_open(self) -> None:
        for script in SCRIPTS:
            subprocess.run([sys.executable, str(script)], check=True)

        data = json.loads(OUT.read_text(encoding="utf-8"))
        self.assertEqual(data["schema_version"], "p2079_s1029_v1")
        self.assertEqual(data["status"], "OPEN_PARTIAL_PROGRESS_WITH_TRACE")
        self.assertEqual(
            data["result_kind"],
            "PASS_BRIDGE_SENSITIVE_IDENTIFIABILITY_AUDIT_WITH_TRACE__C3_STILL_OPEN",
        )

        rows = data["identifiability_results"]["rows"]
        self.assertGreaterEqual(len(rows), 1)
        for row in rows:
            self.assertIn(row["decision"], {"RETAIN_IDENTIFIABLE", "REJECT_DEGENERATE"})
            self.assertGreaterEqual(row["near_optimal_c_span"], 0.0)

        checks = data["gatekeeper_checks"]
        self.assertTrue(checks["preconditions_ready"])
        self.assertTrue(checks["rows_nonempty"])
        self.assertTrue(checks["decision_domain_ok"])
        self.assertFalse(checks["c3_theorem_proven"])
        self.assertTrue(checks["no_background_globalization_claimed"])
        self.assertTrue(checks["no_tensor_component_claimed"])
        self.assertTrue(checks["no_toe_closure_claimed"])


if __name__ == "__main__":
    unittest.main()
