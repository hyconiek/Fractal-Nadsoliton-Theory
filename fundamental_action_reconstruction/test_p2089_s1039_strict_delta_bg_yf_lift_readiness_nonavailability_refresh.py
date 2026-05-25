from __future__ import annotations

import json
import subprocess
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent
SCRIPTS = [
    ROOT / "p2073_s1023_strict_full_ltotal_eom_derivation_scaffold_audit.py",
    ROOT / "p2086_s1036_strict_full_ltotal_eom_termwise_execution_audit.py",
    ROOT / "p2087_s1037_strict_full_ltotal_eom_normal_form_extraction_audit.py",
    ROOT / "p2088_s1038_strict_full_ltotal_eom_theorem_readiness_gap_audit.py",
    ROOT / "p2089_s1039_strict_delta_bg_yf_lift_readiness_nonavailability_refresh.py",
]
OUT = ROOT / "generated" / "p2089_s1039_strict_delta_bg_yf_lift_readiness_nonavailability_refresh.json"


class TestP2089StrictDeltaBgYfLiftReadinessRefresh(unittest.TestCase):
    def test_p2089_exports_next_honest_theorem_object_decision(self) -> None:
        for script in SCRIPTS:
            subprocess.run([sys.executable, str(script)], check=True)

        data = json.loads(OUT.read_text(encoding="utf-8"))
        self.assertEqual(data["schema_version"], "p2089_s1039_v1")
        self.assertEqual(data["status"], "OPEN_PARTIAL_PROGRESS_WITH_TRACE")
        self.assertEqual(
            data["result_kind"],
            "PASS_STRICT_DELTA_BG_YF_LIFT_READINESS_NONAVAILABILITY_REFRESH_WITH_NEXT_STEP__C3_STILL_OPEN",
        )

        checks = data["gatekeeper_checks"]
        self.assertTrue(checks["preconditions_ready"])
        self.assertTrue(checks["r1_local_proxy_pass"])
        self.assertTrue(checks["r2_projection_theorem_open"])
        self.assertTrue(checks["recommended_step_exported"])
        self.assertFalse(checks["c3_theorem_proven"])
        self.assertTrue(checks["no_background_globalization_claimed"])
        self.assertTrue(checks["no_tensor_component_claimed"])
        self.assertTrue(checks["no_toe_closure_claimed"])


if __name__ == "__main__":
    unittest.main()
