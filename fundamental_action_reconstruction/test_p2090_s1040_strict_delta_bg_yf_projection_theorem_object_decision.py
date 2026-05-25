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
    ROOT / "p2090_s1040_strict_delta_bg_yf_projection_theorem_object_decision.py",
]
OUT = ROOT / "generated" / "p2090_s1040_strict_delta_bg_yf_projection_theorem_object_decision.json"


class TestP2090StrictDeltaBgYfProjectionTheoremObjectDecision(unittest.TestCase):
    def test_p2090_exports_current_state_nonavailability_decision(self) -> None:
        for script in SCRIPTS:
            subprocess.run([sys.executable, str(script)], check=True)

        data = json.loads(OUT.read_text(encoding="utf-8"))
        self.assertEqual(data["schema_version"], "p2090_s1040_v1")
        self.assertEqual(data["status"], "OPEN_PARTIAL_PROGRESS_WITH_TRACE")
        self.assertEqual(
            data["result_kind"],
            "PASS_STRICT_DELTA_BG_YF_PROJECTION_THEOREM_OBJECT_DECISION__CURRENT_EXPORT_NONAVAILABILITY__C3_STILL_OPEN",
        )

        checks = data["gatekeeper_checks"]
        self.assertTrue(checks["preconditions_ready"])
        self.assertTrue(checks["decision_mode_is_nonavailability"])
        self.assertTrue(checks["missing_rows_registered"])
        self.assertFalse(checks["c3_theorem_proven"])
        self.assertTrue(checks["no_background_globalization_claimed"])
        self.assertTrue(checks["no_tensor_component_claimed"])
        self.assertTrue(checks["no_toe_closure_claimed"])


if __name__ == "__main__":
    unittest.main()
