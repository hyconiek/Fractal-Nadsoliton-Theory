from __future__ import annotations

import json
import subprocess
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent
SCRIPTS = [
    ROOT / "p2075_s1025_strict_missing_characteristic_candidate_screening_audit.py",
    ROOT / "p2076_s1026_strict_top_candidate_numeric_fit_audit.py",
]
OUT = ROOT / "generated" / "p2076_s1026_strict_top_candidate_numeric_fit_audit.json"


class TestP2076TopCandidateNumericFitAudit(unittest.TestCase):
    def test_p2076_exports_numeric_fit_and_keeps_c3_open(self) -> None:
        for script in SCRIPTS:
            subprocess.run([sys.executable, str(script)], check=True)

        data = json.loads(OUT.read_text(encoding="utf-8"))
        self.assertEqual(data["schema_version"], "p2076_s1026_v1")
        self.assertEqual(data["status"], "OPEN_PARTIAL_PROGRESS_WITH_TRACE")
        self.assertEqual(
            data["result_kind"],
            "PASS_TOP_CANDIDATE_NUMERIC_FIT_AUDIT_WITH_TRACE__C3_STILL_OPEN",
        )

        results = data["fit_results"]["rows"]
        self.assertGreaterEqual(len(results), 1)
        for row in results:
            self.assertIn("absolute_improvement", row)
            self.assertIn(row["decision"], {"RETAIN", "REJECT"})

        checks = data["gatekeeper_checks"]
        self.assertTrue(checks["preconditions_ready"])
        self.assertTrue(checks["fit_rows_nonempty"])
        self.assertTrue(checks["improvement_computed"])
        self.assertTrue(checks["has_retain_or_reject_decisions"])
        self.assertFalse(checks["c3_theorem_proven"])
        self.assertTrue(checks["no_background_globalization_claimed"])
        self.assertTrue(checks["no_tensor_component_claimed"])
        self.assertTrue(checks["no_toe_closure_claimed"])


if __name__ == "__main__":
    unittest.main()
