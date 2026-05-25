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
]
OUT = ROOT / "generated" / "p2086_s1036_strict_full_ltotal_eom_termwise_execution_audit.json"


class TestP2086FullLTotalEOMTermwiseExecutionAudit(unittest.TestCase):
    def test_p2086_exports_termwise_execution_and_keeps_c3_open(self) -> None:
        for script in SCRIPTS:
            subprocess.run([sys.executable, str(script)], check=True)

        data = json.loads(OUT.read_text(encoding="utf-8"))
        self.assertEqual(data["schema_version"], "p2086_s1036_v1")
        self.assertEqual(data["status"], "OPEN_PARTIAL_PROGRESS_WITH_TRACE")
        self.assertEqual(
            data["result_kind"],
            "PASS_STRICT_FULL_LTOTAL_EOM_TERMWISE_EXECUTION_AUDIT_WITH_TRACE__C3_STILL_OPEN",
        )

        checks = data["gatekeeper_checks"]
        self.assertTrue(checks["preconditions_ready"])
        self.assertTrue(checks["all_symbolic_residual_zero"])
        self.assertTrue(checks["all_numeric_probe_residual_zero"])
        self.assertTrue(checks["term_count_expected"])
        self.assertFalse(checks["c3_theorem_proven"])
        self.assertTrue(checks["no_background_globalization_claimed"])
        self.assertTrue(checks["no_tensor_component_claimed"])
        self.assertTrue(checks["no_toe_closure_claimed"])


if __name__ == "__main__":
    unittest.main()
