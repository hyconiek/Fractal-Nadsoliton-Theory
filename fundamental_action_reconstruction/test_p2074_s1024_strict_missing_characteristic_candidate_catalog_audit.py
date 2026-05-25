from __future__ import annotations

import json
import subprocess
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent
SCRIPTS = [
    ROOT / "p2073_s1023_strict_full_ltotal_eom_derivation_scaffold_audit.py",
    ROOT / "p2074_s1024_strict_missing_characteristic_candidate_catalog_audit.py",
]
OUT = ROOT / "generated" / "p2074_s1024_strict_missing_characteristic_candidate_catalog_audit.json"


class TestP2074MissingCharacteristicCandidateCatalogAudit(unittest.TestCase):
    def test_p2074_exports_candidate_catalog_and_keeps_c3_open(self) -> None:
        for script in SCRIPTS:
            subprocess.run([sys.executable, str(script)], check=True)

        data = json.loads(OUT.read_text(encoding="utf-8"))
        self.assertEqual(data["schema_version"], "p2074_s1024_v1")
        self.assertEqual(data["status"], "OPEN_PARTIAL_PROGRESS_WITH_TRACE")
        self.assertEqual(
            data["result_kind"],
            "PASS_MISSING_CHARACTERISTIC_CANDIDATE_CATALOG_AUDIT_WITH_TRACE__C3_STILL_OPEN",
        )

        cat = data["candidate_catalog"]
        self.assertGreaterEqual(cat["candidate_count"], 1)
        self.assertEqual(cat["candidate_count"], len(cat["rows"]))
        for row in cat["rows"]:
            self.assertIn("delta_eom", row)
            self.assertIn("selector_qw2191_screening", row)

        checks = data["gatekeeper_checks"]
        self.assertTrue(checks["preconditions_ready"])
        self.assertTrue(checks["candidate_rows_nonempty"])
        self.assertTrue(checks["delta_eom_exported"])
        self.assertTrue(checks["qw2191_proxy_screening_exported"])
        self.assertFalse(checks["c3_theorem_proven"])
        self.assertTrue(checks["no_background_globalization_claimed"])
        self.assertTrue(checks["no_tensor_component_claimed"])
        self.assertTrue(checks["no_toe_closure_claimed"])


if __name__ == "__main__":
    unittest.main()
