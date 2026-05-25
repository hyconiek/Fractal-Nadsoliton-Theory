from __future__ import annotations

import json
import subprocess
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent
SCRIPTS = [
    ROOT / "p2076_s1026_strict_top_candidate_numeric_fit_audit.py",
    ROOT / "p2077_s1027_strict_constrained_coefficient_fit_audit.py",
]
OUT = ROOT / "generated" / "p2077_s1027_strict_constrained_coefficient_fit_audit.json"


class TestP2077ConstrainedCoefficientFitAudit(unittest.TestCase):
    def test_p2077_exports_constrained_fit_and_keeps_c3_open(self) -> None:
        for script in SCRIPTS:
            subprocess.run([sys.executable, str(script)], check=True)

        data = json.loads(OUT.read_text(encoding="utf-8"))
        self.assertEqual(data["schema_version"], "p2077_s1027_v1")
        self.assertEqual(data["status"], "OPEN_PARTIAL_PROGRESS_WITH_TRACE")
        self.assertEqual(
            data["result_kind"],
            "PASS_CONSTRAINED_COEFFICIENT_FIT_AUDIT_WITH_TRACE__C3_STILL_OPEN",
        )

        rows = data["fit_results"]["rows"]
        self.assertGreaterEqual(len(rows), 1)
        for row in rows:
            self.assertLessEqual(abs(float(row["best_coefficient"])), 0.2)
            self.assertEqual(row["qw2191_proxy_risk"], "LOW_PROXY")
            self.assertIn(row["decision"], {"RETAIN", "REJECT"})

        checks = data["gatekeeper_checks"]
        self.assertTrue(checks["preconditions_ready"])
        self.assertTrue(checks["rows_nonempty"])
        self.assertTrue(checks["all_coefficients_within_bound"])
        self.assertTrue(checks["all_qw2191_proxy_low"])
        self.assertFalse(checks["c3_theorem_proven"])
        self.assertTrue(checks["no_background_globalization_claimed"])
        self.assertTrue(checks["no_tensor_component_claimed"])
        self.assertTrue(checks["no_toe_closure_claimed"])


if __name__ == "__main__":
    unittest.main()
