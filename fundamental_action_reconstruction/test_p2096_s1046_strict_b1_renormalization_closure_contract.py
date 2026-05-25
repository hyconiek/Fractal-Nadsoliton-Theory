from __future__ import annotations

import json
import subprocess
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent
SCRIPTS = [
    ROOT / "p1950_s900_strict_renormalization_exact_integration.py",
    ROOT / "p2094_s1044_strict_b1_quotient_renormalization_rank_repair.py",
    ROOT / "p2095_s1045_strict_b1_gb_derived_channel_certificate.py",
    ROOT / "p2096_s1046_strict_b1_renormalization_closure_contract.py",
]
OUT = ROOT / "generated" / "p2096_s1046_strict_b1_renormalization_closure_contract.json"


class TestP2096StrictB1RenormalizationClosureContract(unittest.TestCase):
    def test_p2096_exports_quotient_scope_contract(self) -> None:
        for script in SCRIPTS:
            subprocess.run([sys.executable, str(script)], check=True)

        data = json.loads(OUT.read_text(encoding="utf-8"))
        self.assertEqual(data["schema_version"], "p2096_s1046_v1")
        self.assertEqual(data["status"], "OPEN_PARTIAL_PROGRESS_WITH_TRACE")
        self.assertEqual(
            data["result_kind"],
            "PASS_STRICT_B1_RENORMALIZATION_CLOSURE_CONTRACT_WITH_TRACE__QUOTIENT_SCOPE_ONLY",
        )

        checks = data["gatekeeper_checks"]
        self.assertTrue(checks["preconditions_ready"])
        self.assertTrue(checks["quotient_scope_contract_exported"])
        self.assertFalse(checks["full_4channel_independence_proven"])
        self.assertFalse(checks["global_background_family_closure_proven"])
        self.assertFalse(checks["c3_theorem_proven"])


if __name__ == "__main__":
    unittest.main()
