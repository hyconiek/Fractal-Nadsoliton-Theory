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
    ROOT / "p2128_s1078_strict_cmp2_block_bootstrap_dependence_aware_uncertainty_inflation.py",
]
OUT = ROOT / "generated" / "p2128_s1078_strict_cmp2_block_bootstrap_dependence_aware_uncertainty_inflation.json"


class TestP2128StrictCmp2BlockBootstrapDependenceAwareUncertaintyInflation(unittest.TestCase):
    def test_p2128_exports_dependence_aware_inflation_object(self) -> None:
        for script in SCRIPTS:
            subprocess.run([sys.executable, str(script)], check=True)

        data = json.loads(OUT.read_text(encoding="utf-8"))
        self.assertEqual(data["schema_version"], "p2128_s1078_v1")
        self.assertEqual(data["status"], "OPEN_PARTIAL_PROGRESS_WITH_TRACE")
        self.assertEqual(
            data["result_kind"],
            "PASS_STRICT_CMP2_BLOCK_BOOTSTRAP_DEPENDENCE_AWARE_UNCERTAINTY_INFLATION_WITH_TRACE",
        )

        checks = data["gatekeeper_checks"]
        self.assertTrue(checks["preconditions_ready"])
        self.assertTrue(checks["block_bootstrap_executed"])
        self.assertTrue(checks["inflation_object_exported"])
        self.assertFalse(checks["full_d3_covariance_transport_proven"])
        self.assertFalse(checks["full_cutkosky_closure_proven"])
        self.assertFalse(checks["c3_theorem_proven"])
        self.assertTrue(checks["no_toe_closure_claimed"])


if __name__ == "__main__":
    unittest.main()
