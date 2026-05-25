from __future__ import annotations

import json
import subprocess
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent
SCRIPTS = [
    ROOT / "p2069_s1019_strict_same_scheme_tau_remediation_scenario_audit.py",
    ROOT / "p2070_s1020_strict_same_scheme_tau_grid_extension_feasibility_audit.py",
]
OUT = ROOT / "generated" / "p2070_s1020_strict_same_scheme_tau_grid_extension_feasibility_audit.json"


class TestP2070TauGridExtensionFeasibilityAudit(unittest.TestCase):
    def test_p2070_exports_extended_grid_and_keeps_c3_open(self) -> None:
        for script in SCRIPTS:
            subprocess.run([sys.executable, str(script)], check=True)

        data = json.loads(OUT.read_text(encoding="utf-8"))
        self.assertEqual(data["schema_version"], "p2070_s1020_v1")
        self.assertEqual(data["status"], "OPEN_PARTIAL_PROGRESS_WITH_TRACE")
        self.assertEqual(
            data["result_kind"],
            "PASS_TAU_GRID_EXTENSION_FEASIBILITY_AUDIT_WITH_TRACE__C3_STILL_OPEN",
        )

        tau = data["tau_grid_extension_audit"]
        self.assertGreaterEqual(len(tau["extension_rows"]), 1)
        self.assertIsNotNone(tau["maximal_feasible_tau_in_extension"])
        self.assertGreaterEqual(tau["required_tau_reduction_factor"], 0.0)
        self.assertLessEqual(tau["required_tau_reduction_factor"], 1.0)

        checks = data["gatekeeper_checks"]
        self.assertTrue(checks["preconditions_ready"])
        self.assertTrue(checks["extension_rows_nonempty"])
        self.assertTrue(checks["maximal_feasible_exists"])
        self.assertTrue(checks["required_tau_reduction_factor_in_0_1"])
        self.assertTrue(checks["required_tau_reduction_percent_nonnegative"])
        self.assertFalse(checks["c3_theorem_proven"])
        self.assertTrue(checks["no_background_globalization_claimed"])
        self.assertTrue(checks["no_tensor_component_claimed"])
        self.assertTrue(checks["no_toe_closure_claimed"])


if __name__ == "__main__":
    unittest.main()
