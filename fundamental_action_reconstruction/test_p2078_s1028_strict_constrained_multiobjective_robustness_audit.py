from __future__ import annotations

import json
import subprocess
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent
SCRIPTS = [
    ROOT / "p2077_s1027_strict_constrained_coefficient_fit_audit.py",
    ROOT / "p2078_s1028_strict_constrained_multiobjective_robustness_audit.py",
]
OUT = ROOT / "generated" / "p2078_s1028_strict_constrained_multiobjective_robustness_audit.json"


class TestP2078ConstrainedMultiobjectiveRobustnessAudit(unittest.TestCase):
    def test_p2078_exports_robustness_and_keeps_c3_open(self) -> None:
        for script in SCRIPTS:
            subprocess.run([sys.executable, str(script)], check=True)

        data = json.loads(OUT.read_text(encoding="utf-8"))
        self.assertEqual(data["schema_version"], "p2078_s1028_v1")
        self.assertEqual(data["status"], "OPEN_PARTIAL_PROGRESS_WITH_TRACE")
        self.assertEqual(
            data["result_kind"],
            "PASS_CONSTRAINED_MULTIOBJECTIVE_ROBUSTNESS_AUDIT_WITH_TRACE__C3_STILL_OPEN",
        )

        rows = data["robustness_results"]["rows"]
        self.assertGreaterEqual(len(rows), 1)
        for row in rows:
            self.assertIn(row["decision"], {"RETAIN_ROBUST", "REJECT_FRAGILE"})
            self.assertGreaterEqual(row["robustness_ratio"], 0.0)
            self.assertLessEqual(row["robustness_ratio"], 1.0)

        checks = data["gatekeeper_checks"]
        self.assertTrue(checks["preconditions_ready"])
        self.assertTrue(checks["candidate_rows_nonempty"])
        self.assertTrue(checks["has_robust_or_fragile_decisions"])
        self.assertFalse(checks["c3_theorem_proven"])
        self.assertTrue(checks["no_background_globalization_claimed"])
        self.assertTrue(checks["no_tensor_component_claimed"])
        self.assertTrue(checks["no_toe_closure_claimed"])


if __name__ == "__main__":
    unittest.main()
