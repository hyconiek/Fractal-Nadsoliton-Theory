from __future__ import annotations

import json
import subprocess
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent
SCRIPTS = [
    ROOT / "p2126_s1076_strict_cmp2_backend_evidence_weighted_posterior_calibration_audit.py",
    ROOT / "p2127_s1077_strict_cmp2_bootstrap_backend_evidence_stresstest.py",
    ROOT / "p2129_s1079_strict_cmp2_block_definition_sensitivity_sweep.py",
    ROOT / "p2130_s1080_strict_cmp2_sample_support_expansion_planner.py",
]
OUT = ROOT / "generated" / "p2130_s1080_strict_cmp2_sample_support_expansion_planner.json"


class TestP2130StrictCmp2SampleSupportExpansionPlanner(unittest.TestCase):
    def test_p2130_exports_support_expansion_plan(self) -> None:
        for script in SCRIPTS:
            subprocess.run([sys.executable, str(script)], check=True)

        data = json.loads(OUT.read_text(encoding="utf-8"))
        self.assertEqual(data["schema_version"], "p2130_s1080_v1")
        self.assertEqual(data["status"], "OPEN_PARTIAL_PROGRESS_WITH_TRACE")
        self.assertEqual(
            data["result_kind"],
            "PASS_STRICT_CMP2_SAMPLE_SUPPORT_EXPANSION_PLANNER_WITH_TRACE",
        )

        checks = data["gatekeeper_checks"]
        self.assertTrue(checks["preconditions_ready"])
        self.assertTrue(checks["root_cause_diagnosed"])
        self.assertTrue(checks["support_expansion_plan_exported"])
        self.assertFalse(checks["full_d3_covariance_transport_proven"])
        self.assertFalse(checks["full_cutkosky_closure_proven"])
        self.assertFalse(checks["c3_theorem_proven"])
        self.assertTrue(checks["no_toe_closure_claimed"])

        plan = data["sample_support_expansion_planner"]["projected_ci95_widths_by_support_multiplier"]
        self.assertIn("x2", plan)
        self.assertIn("x3", plan)
        self.assertIn("x5", plan)


if __name__ == "__main__":
    unittest.main()
