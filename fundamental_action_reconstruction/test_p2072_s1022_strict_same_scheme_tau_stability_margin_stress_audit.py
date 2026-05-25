from __future__ import annotations

import json
import subprocess
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent
SCRIPTS = [
    ROOT / "p2071_s1021_strict_same_scheme_tau_boundary_refinement_audit.py",
    ROOT / "p2072_s1022_strict_same_scheme_tau_stability_margin_stress_audit.py",
]
OUT = ROOT / "generated" / "p2072_s1022_strict_same_scheme_tau_stability_margin_stress_audit.json"


class TestP2072TauStabilityMarginStressAudit(unittest.TestCase):
    def test_p2072_exports_stress_margin_and_keeps_c3_open(self) -> None:
        for script in SCRIPTS:
            subprocess.run([sys.executable, str(script)], check=True)

        data = json.loads(OUT.read_text(encoding="utf-8"))
        self.assertEqual(data["schema_version"], "p2072_s1022_v1")
        self.assertEqual(data["status"], "OPEN_PARTIAL_PROGRESS_WITH_TRACE")
        self.assertEqual(
            data["result_kind"],
            "PASS_TAU_STABILITY_MARGIN_STRESS_AUDIT_WITH_TRACE__C3_STILL_OPEN",
        )

        m = data["tau_stability_margin"]
        self.assertGreaterEqual(len(m["stress_rows"]), 1)
        self.assertIn(m["boundary_classification"], {"SHARP_STABLE_BOUNDARY", "FRAGILE_OR_UNRESOLVED_BOUNDARY"})

        checks = data["gatekeeper_checks"]
        self.assertTrue(checks["preconditions_ready"])
        self.assertTrue(checks["stress_rows_nonempty"])
        self.assertTrue(checks["contains_requested_deltas"])
        self.assertFalse(checks["c3_theorem_proven"])
        self.assertTrue(checks["no_background_globalization_claimed"])
        self.assertTrue(checks["no_tensor_component_claimed"])
        self.assertTrue(checks["no_toe_closure_claimed"])


if __name__ == "__main__":
    unittest.main()
