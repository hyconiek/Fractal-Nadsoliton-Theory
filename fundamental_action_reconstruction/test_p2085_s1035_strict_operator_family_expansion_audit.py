from __future__ import annotations

import json
import subprocess
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent
SCRIPTS = [
    ROOT / "p2084_s1034_strict_adaptive_channel_weight_audit.py",
    ROOT / "p2085_s1035_strict_operator_family_expansion_audit.py",
]
OUT = ROOT / "generated" / "p2085_s1035_strict_operator_family_expansion_audit.json"


class TestP2085OperatorFamilyExpansionAudit(unittest.TestCase):
    def test_p2085_exports_expansion_audit_and_keeps_c3_open(self) -> None:
        for script in SCRIPTS:
            subprocess.run([sys.executable, str(script)], check=True)

        data = json.loads(OUT.read_text(encoding="utf-8"))
        self.assertEqual(data["schema_version"], "p2085_s1035_v1")
        self.assertEqual(data["status"], "OPEN_PARTIAL_PROGRESS_WITH_TRACE")
        self.assertEqual(
            data["result_kind"],
            "PASS_OPERATOR_FAMILY_EXPANSION_AUDIT_WITH_TRACE__C3_STILL_OPEN",
        )

        res = data["expansion_results"]
        self.assertEqual(len(res["family_results"]), 2)
        self.assertGreaterEqual(res["robust_family_count"], 0)

        checks = data["gatekeeper_checks"]
        self.assertTrue(checks["preconditions_ready"])
        self.assertTrue(checks["family_count_expected"])
        self.assertTrue(checks["all_coefficients_within_bound"])
        self.assertTrue(checks["decision_domain_ok"])
        self.assertFalse(checks["c3_theorem_proven"])
        self.assertTrue(checks["no_background_globalization_claimed"])
        self.assertTrue(checks["no_tensor_component_claimed"])
        self.assertTrue(checks["no_toe_closure_claimed"])


if __name__ == "__main__":
    unittest.main()
