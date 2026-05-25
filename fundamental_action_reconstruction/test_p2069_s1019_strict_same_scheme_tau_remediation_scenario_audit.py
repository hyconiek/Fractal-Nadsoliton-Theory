from __future__ import annotations

import json
import subprocess
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent
SCRIPTS = [
    ROOT / "p2068_s1018_strict_same_scheme_transportability_improvement_lever_audit.py",
    ROOT / "p2069_s1019_strict_same_scheme_tau_remediation_scenario_audit.py",
]
OUT = ROOT / "generated" / "p2069_s1019_strict_same_scheme_tau_remediation_scenario_audit.json"


class TestP2069TauRemediationScenarioAudit(unittest.TestCase):
    def test_p2069_exports_tau_scenarios_and_keeps_c3_open(self) -> None:
        for script in SCRIPTS:
            subprocess.run([sys.executable, str(script)], check=True)

        data = json.loads(OUT.read_text(encoding="utf-8"))
        self.assertEqual(data["schema_version"], "p2069_s1019_v1")
        self.assertEqual(data["status"], "OPEN_PARTIAL_PROGRESS_WITH_TRACE")
        self.assertEqual(
            data["result_kind"],
            "PASS_TAU_REMEDIATION_SCENARIO_AUDIT_WITH_TRACE__C3_STILL_OPEN",
        )

        tau = data["tau_remediation"]
        self.assertGreaterEqual(tau["scenario_count"], 1)
        self.assertEqual(tau["scenario_count"], len(tau["scenario_rows"]))

        checks = data["gatekeeper_checks"]
        self.assertTrue(checks["preconditions_ready"])
        self.assertTrue(checks["scenario_rows_nonempty"])
        self.assertTrue(checks["required_tau_reduction_factor_in_0_1"])
        self.assertTrue(checks["required_tau_reduction_percent_nonnegative"])
        self.assertFalse(checks["c3_theorem_proven"])
        self.assertTrue(checks["no_background_globalization_claimed"])
        self.assertTrue(checks["no_tensor_component_claimed"])
        self.assertTrue(checks["no_toe_closure_claimed"])


if __name__ == "__main__":
    unittest.main()
