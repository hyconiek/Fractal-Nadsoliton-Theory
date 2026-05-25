from __future__ import annotations

import json
import subprocess
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent
SCRIPTS = [
    ROOT / "p2072_s1022_strict_same_scheme_tau_stability_margin_stress_audit.py",
    ROOT / "p2073_s1023_strict_full_ltotal_eom_derivation_scaffold_audit.py",
]
OUT = ROOT / "generated" / "p2073_s1023_strict_full_ltotal_eom_derivation_scaffold_audit.json"


class TestP2073StrictFullLTotalEOMDerivationScaffoldAudit(unittest.TestCase):
    def test_p2073_exports_scaffold_and_keeps_c3_open(self) -> None:
        for script in SCRIPTS:
            subprocess.run([sys.executable, str(script)], check=True)

        data = json.loads(OUT.read_text(encoding="utf-8"))
        self.assertEqual(data["schema_version"], "p2073_s1023_v1")
        self.assertEqual(data["status"], "OPEN_PARTIAL_PROGRESS_WITH_TRACE")
        self.assertEqual(
            data["result_kind"],
            "PASS_STRICT_FULL_LTOTAL_EOM_DERIVATION_SCAFFOLD_AUDIT_WITH_TRACE__C3_STILL_OPEN",
        )

        kctx = data["kernel_split_context"]
        self.assertEqual(kctx["bridge_status"], "OPEN_NO_RIGOROUS_IDENTIFICATION_EXPORT")
        self.assertEqual(kctx["missing_characteristic_hypothesis"], "OPEN_CANDIDATE_REQUIRED")

        scaffold = data["eom_scaffold"]
        self.assertGreaterEqual(len(scaffold["fields"]), 1)
        self.assertGreaterEqual(len(scaffold["variation_map"]), 1)

        checks = data["gatekeeper_checks"]
        self.assertTrue(checks["preconditions_ready"])
        self.assertTrue(checks["scaffold_fields_nonempty"])
        self.assertTrue(checks["variation_map_nonempty"])
        self.assertTrue(checks["symbolic_trace_scaffold_only"])
        self.assertFalse(checks["c3_theorem_proven"])
        self.assertTrue(checks["no_background_globalization_claimed"])
        self.assertTrue(checks["no_tensor_component_claimed"])
        self.assertTrue(checks["no_toe_closure_claimed"])


if __name__ == "__main__":
    unittest.main()
