from __future__ import annotations

import json
import subprocess
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent
SCRIPTS = [
    ROOT / "p2083_s1033_strict_three_feature_targeted_composite_operator_audit.py",
    ROOT / "p2084_s1034_strict_adaptive_channel_weight_audit.py",
]
OUT = ROOT / "generated" / "p2084_s1034_strict_adaptive_channel_weight_audit.json"


class TestP2084AdaptiveChannelWeightAudit(unittest.TestCase):
    def test_p2084_exports_adaptive_weight_audit_and_keeps_c3_open(self) -> None:
        for script in SCRIPTS:
            subprocess.run([sys.executable, str(script)], check=True)

        data = json.loads(OUT.read_text(encoding="utf-8"))
        self.assertEqual(data["schema_version"], "p2084_s1034_v1")
        self.assertEqual(data["status"], "OPEN_PARTIAL_PROGRESS_WITH_TRACE")
        self.assertEqual(
            data["result_kind"],
            "PASS_ADAPTIVE_CHANNEL_WEIGHT_AUDIT_WITH_TRACE__C3_STILL_OPEN",
        )

        res = data["adaptive_results"]
        self.assertIn(res["decision"], {"RETAIN_ADAPTIVE_REGION", "REJECT_NO_STABLE_REGION"})
        self.assertGreaterEqual(res["feasible_region_size"], 0)

        checks = data["gatekeeper_checks"]
        self.assertTrue(checks["preconditions_ready"])
        self.assertTrue(checks["decision_domain_ok"])
        self.assertTrue(checks["global_best_coefficients_within_bound"])
        self.assertFalse(checks["c3_theorem_proven"])
        self.assertTrue(checks["no_background_globalization_claimed"])
        self.assertTrue(checks["no_tensor_component_claimed"])
        self.assertTrue(checks["no_toe_closure_claimed"])


if __name__ == "__main__":
    unittest.main()
