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
]
OUT = ROOT / "generated" / "p2129_s1079_strict_cmp2_block_definition_sensitivity_sweep.json"


class TestP2129StrictCmp2BlockDefinitionSensitivitySweep(unittest.TestCase):
    def test_p2129_exports_block_sensitivity_sweep(self) -> None:
        for script in SCRIPTS:
            subprocess.run([sys.executable, str(script)], check=True)

        data = json.loads(OUT.read_text(encoding="utf-8"))
        self.assertEqual(data["schema_version"], "p2129_s1079_v1")
        self.assertEqual(data["status"], "OPEN_PARTIAL_PROGRESS_WITH_TRACE")
        self.assertEqual(
            data["result_kind"],
            "PASS_STRICT_CMP2_BLOCK_DEFINITION_SENSITIVITY_SWEEP_WITH_TRACE",
        )

        checks = data["gatekeeper_checks"]
        self.assertTrue(checks["preconditions_ready"])
        self.assertTrue(checks["sweep_executed"])
        self.assertTrue(checks["inflation_stability_compared"])
        self.assertFalse(checks["full_d3_covariance_transport_proven"])
        self.assertFalse(checks["full_cutkosky_closure_proven"])
        self.assertFalse(checks["c3_theorem_proven"])
        self.assertTrue(checks["no_toe_closure_claimed"])

        variants = data["block_definition_sensitivity_sweep"]["variants"]
        self.assertIn("backend_s_unique", variants)
        self.assertIn("backend_s_quantile_q2", variants)
        self.assertIn("backend_s_quantile_q4", variants)
        self.assertIn("bin_adjacency_k2", variants)


if __name__ == "__main__":
    unittest.main()
