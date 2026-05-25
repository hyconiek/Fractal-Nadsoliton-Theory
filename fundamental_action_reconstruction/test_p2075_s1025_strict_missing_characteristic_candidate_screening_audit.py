from __future__ import annotations

import json
import subprocess
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent
SCRIPTS = [
    ROOT / "p2074_s1024_strict_missing_characteristic_candidate_catalog_audit.py",
    ROOT / "p2075_s1025_strict_missing_characteristic_candidate_screening_audit.py",
]
OUT = ROOT / "generated" / "p2075_s1025_strict_missing_characteristic_candidate_screening_audit.json"


class TestP2075MissingCharacteristicCandidateScreeningAudit(unittest.TestCase):
    def test_p2075_exports_screening_and_keeps_c3_open(self) -> None:
        for script in SCRIPTS:
            subprocess.run([sys.executable, str(script)], check=True)

        data = json.loads(OUT.read_text(encoding="utf-8"))
        self.assertEqual(data["schema_version"], "p2075_s1025_v1")
        self.assertEqual(data["status"], "OPEN_PARTIAL_PROGRESS_WITH_TRACE")
        self.assertEqual(
            data["result_kind"],
            "PASS_MISSING_CHARACTERISTIC_CANDIDATE_SCREENING_AUDIT_WITH_TRACE__C3_STILL_OPEN",
        )

        table = data["screening_table"]
        self.assertGreaterEqual(table["candidate_count"], 1)
        self.assertEqual(table["candidate_count"], len(table["rows"]))

        checks = data["gatekeeper_checks"]
        self.assertTrue(checks["preconditions_ready"])
        self.assertTrue(checks["screening_rows_nonempty"])
        self.assertTrue(checks["priority_ranking_monotone"])
        self.assertFalse(checks["c3_theorem_proven"])
        self.assertTrue(checks["no_background_globalization_claimed"])
        self.assertTrue(checks["no_tensor_component_claimed"])
        self.assertTrue(checks["no_toe_closure_claimed"])


if __name__ == "__main__":
    unittest.main()
