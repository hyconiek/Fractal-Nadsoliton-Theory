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
]
OUT = ROOT / "generated" / "p2094_s1044_strict_b1_quotient_renormalization_rank_repair.json"


class TestP2094StrictB1QuotientRenormalizationRankRepair(unittest.TestCase):
    def test_p2094_exports_rank_repair_witness(self) -> None:
        for script in SCRIPTS:
            subprocess.run([sys.executable, str(script)], check=True)

        data = json.loads(OUT.read_text(encoding="utf-8"))
        self.assertEqual(data["schema_version"], "p2094_s1044_v1")
        self.assertEqual(data["status"], "OPEN_PARTIAL_PROGRESS_WITH_TRACE")
        self.assertEqual(
            data["result_kind"],
            "PASS_STRICT_B1_QUOTIENT_RENORMALIZATION_RANK_REPAIR_WITH_TRACE__NO_FULL_GB_INDEPENDENCE",
        )

        checks = data["gatekeeper_checks"]
        self.assertTrue(checks["preconditions_ready"])
        self.assertTrue(checks["quotient_3x3_residual_small"])
        self.assertTrue(checks["gb_row_consistent_as_derived_channel"])
        self.assertFalse(checks["full_4channel_independence_proven"])
        self.assertFalse(checks["c3_theorem_proven"])


if __name__ == "__main__":
    unittest.main()
