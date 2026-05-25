from __future__ import annotations

import json
import subprocess
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent
SCRIPTS = [
    ROOT / "p2070_s1020_strict_same_scheme_tau_grid_extension_feasibility_audit.py",
    ROOT / "p2071_s1021_strict_same_scheme_tau_boundary_refinement_audit.py",
]
OUT = ROOT / "generated" / "p2071_s1021_strict_same_scheme_tau_boundary_refinement_audit.json"


class TestP2071TauBoundaryRefinementAudit(unittest.TestCase):
    def test_p2071_exports_boundary_refinement_and_keeps_c3_open(self) -> None:
        for script in SCRIPTS:
            subprocess.run([sys.executable, str(script)], check=True)

        data = json.loads(OUT.read_text(encoding="utf-8"))
        self.assertEqual(data["schema_version"], "p2071_s1021_v1")
        self.assertEqual(data["status"], "OPEN_PARTIAL_PROGRESS_WITH_TRACE")
        self.assertEqual(data["result_kind"], "PASS_TAU_BOUNDARY_REFINEMENT_AUDIT_WITH_TRACE__C3_STILL_OPEN")

        ref = data["tau_boundary_refinement"]
        self.assertGreaterEqual(len(ref["rows"]), 1)
        self.assertIsNotNone(ref["tight_feasible_tau"])
        self.assertGreaterEqual(ref["required_tau_reduction_factor"], 0.0)
        self.assertLessEqual(ref["required_tau_reduction_factor"], 1.0)

        checks = data["gatekeeper_checks"]
        self.assertTrue(checks["preconditions_ready"])
        self.assertTrue(checks["rows_nonempty"])
        self.assertTrue(checks["tight_feasible_exists"])
        self.assertTrue(checks["required_tau_reduction_factor_in_0_1"])
        self.assertFalse(checks["c3_theorem_proven"])
        self.assertTrue(checks["no_background_globalization_claimed"])
        self.assertTrue(checks["no_tensor_component_claimed"])
        self.assertTrue(checks["no_toe_closure_claimed"])


if __name__ == "__main__":
    unittest.main()
