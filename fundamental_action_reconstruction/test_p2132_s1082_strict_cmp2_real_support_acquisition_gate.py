from __future__ import annotations

import json
import subprocess
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent
SCRIPTS = [
    ROOT / "p2126_s1076_strict_cmp2_backend_evidence_weighted_posterior_calibration_audit.py",
    ROOT / "p2132_s1082_strict_cmp2_real_support_acquisition_gate.py",
]
OUT = ROOT / "generated" / "p2132_s1082_strict_cmp2_real_support_acquisition_gate.json"


class TestP2132StrictCmp2RealSupportAcquisitionGate(unittest.TestCase):
    def test_p2132_exports_readiness_gate(self) -> None:
        for script in SCRIPTS:
            subprocess.run([sys.executable, str(script)], check=True)

        data = json.loads(OUT.read_text(encoding="utf-8"))
        self.assertEqual(data["schema_version"], "p2132_s1082_v1")
        self.assertEqual(data["status"], "OPEN_PARTIAL_PROGRESS_WITH_TRACE")
        self.assertIn(data["result_kind"], {
            "PASS_STRICT_CMP2_REAL_SUPPORT_ACQUISITION_GATE_WITH_TRACE",
            "OPEN_STRICT_CMP2_REAL_SUPPORT_ACQUISITION_GATE_BLOCKED",
        })

        checks = data["gatekeeper_checks"]
        self.assertTrue(checks["base_rows_present"])
        self.assertFalse(checks["full_d3_covariance_transport_proven"])
        self.assertFalse(checks["full_cutkosky_closure_proven"])
        self.assertFalse(checks["c3_theorem_proven"])
        self.assertTrue(checks["no_toe_closure_claimed"])


if __name__ == "__main__":
    unittest.main()
