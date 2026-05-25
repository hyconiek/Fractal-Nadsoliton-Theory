from __future__ import annotations

import json
import subprocess
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent
SCRIPTS = [
    ROOT / "p2080_s1030_strict_profile_shape_bridge_audit.py",
    ROOT / "p2081_s1031_strict_feature_weighted_candidate_redesign_audit.py",
]
OUT = ROOT / "generated" / "p2081_s1031_strict_feature_weighted_candidate_redesign_audit.json"


class TestP2081FeatureWeightedCandidateRedesignAudit(unittest.TestCase):
    def test_p2081_exports_redesign_and_keeps_c3_open(self) -> None:
        for script in SCRIPTS:
            subprocess.run([sys.executable, str(script)], check=True)

        data = json.loads(OUT.read_text(encoding="utf-8"))
        self.assertEqual(data["schema_version"], "p2081_s1031_v1")
        self.assertEqual(data["status"], "OPEN_PARTIAL_PROGRESS_WITH_TRACE")
        self.assertEqual(
            data["result_kind"],
            "PASS_FEATURE_WEIGHTED_CANDIDATE_REDESIGN_AUDIT_WITH_TRACE__C3_STILL_OPEN",
        )

        rows = data["redesign_results"]["rows"]
        self.assertGreaterEqual(len(rows), 1)
        for row in rows:
            self.assertLessEqual(abs(float(row["best_coefficient"])), 0.2)
            self.assertIn(row["decision"], {"RETAIN_REDESIGNED", "REJECT_REDESIGNED"})

        checks = data["gatekeeper_checks"]
        self.assertTrue(checks["preconditions_ready"])
        self.assertTrue(checks["rows_nonempty"])
        self.assertTrue(checks["coefficients_within_bound"])
        self.assertTrue(checks["decision_domain_ok"])
        self.assertFalse(checks["c3_theorem_proven"])
        self.assertTrue(checks["no_background_globalization_claimed"])
        self.assertTrue(checks["no_tensor_component_claimed"])
        self.assertTrue(checks["no_toe_closure_claimed"])


if __name__ == "__main__":
    unittest.main()
