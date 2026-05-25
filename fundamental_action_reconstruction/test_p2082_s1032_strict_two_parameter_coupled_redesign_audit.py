from __future__ import annotations

import json
import subprocess
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent
SCRIPTS = [
    ROOT / "p2081_s1031_strict_feature_weighted_candidate_redesign_audit.py",
    ROOT / "p2082_s1032_strict_two_parameter_coupled_redesign_audit.py",
]
OUT = ROOT / "generated" / "p2082_s1032_strict_two_parameter_coupled_redesign_audit.json"


class TestP2082TwoParameterCoupledRedesignAudit(unittest.TestCase):
    def test_p2082_exports_coupled_redesign_and_keeps_c3_open(self) -> None:
        for script in SCRIPTS:
            subprocess.run([sys.executable, str(script)], check=True)

        data = json.loads(OUT.read_text(encoding="utf-8"))
        self.assertEqual(data["schema_version"], "p2082_s1032_v1")
        self.assertEqual(data["status"], "OPEN_PARTIAL_PROGRESS_WITH_TRACE")
        self.assertEqual(
            data["result_kind"],
            "PASS_TWO_PARAMETER_COUPLED_REDESIGN_AUDIT_WITH_TRACE__C3_STILL_OPEN",
        )

        res = data["coupled_results"]
        self.assertLessEqual(abs(float(res["best_c1"])), 0.2)
        self.assertLessEqual(abs(float(res["best_c2"])), 0.2)
        self.assertIn(res["decision"], {"RETAIN_COUPLED", "REJECT_COUPLED"})

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
