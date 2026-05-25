from __future__ import annotations

import json
import subprocess
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent
SCRIPTS = [
    ROOT / "p2082_s1032_strict_two_parameter_coupled_redesign_audit.py",
    ROOT / "p2083_s1033_strict_three_feature_targeted_composite_operator_audit.py",
]
OUT = ROOT / "generated" / "p2083_s1033_strict_three_feature_targeted_composite_operator_audit.json"


class TestP2083ThreeFeatureTargetedCompositeOperatorAudit(unittest.TestCase):
    def test_p2083_exports_composite_audit_and_keeps_c3_open(self) -> None:
        for script in SCRIPTS:
            subprocess.run([sys.executable, str(script)], check=True)

        data = json.loads(OUT.read_text(encoding="utf-8"))
        self.assertEqual(data["schema_version"], "p2083_s1033_v1")
        self.assertEqual(data["status"], "OPEN_PARTIAL_PROGRESS_WITH_TRACE")
        self.assertEqual(
            data["result_kind"],
            "PASS_THREE_FEATURE_TARGETED_COMPOSITE_OPERATOR_AUDIT_WITH_TRACE__C3_STILL_OPEN",
        )

        res = data["composite_results"]
        self.assertLessEqual(abs(float(res["best_c1"])), 0.2)
        self.assertLessEqual(abs(float(res["best_c2"])), 0.2)
        self.assertLessEqual(abs(float(res["best_c3"])), 0.2)
        self.assertIn(res["decision"], {"RETAIN_COMPOSITE", "REJECT_COMPOSITE"})

        checks = data["gatekeeper_checks"]
        self.assertTrue(checks["preconditions_ready"])
        self.assertTrue(checks["coefficients_within_bound"])
        self.assertTrue(checks["decision_domain_ok"])
        self.assertFalse(checks["c3_theorem_proven"])
        self.assertTrue(checks["no_background_globalization_claimed"])
        self.assertTrue(checks["no_tensor_component_claimed"])
        self.assertTrue(checks["no_toe_closure_claimed"])


if __name__ == "__main__":
    unittest.main()
