from __future__ import annotations

import json
import subprocess
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent
P2037 = ROOT / "p2037_s987_strict_task1_same_scheme_finite_part_map_seed.py"
P2038 = ROOT / "p2038_s988_strict_same_scheme_finite_part_map_data_import_audit.py"
P2039 = ROOT / "p2039_s989_strict_same_scheme_finite_part_candidate_uncertainty_bound_probe.py"
P2040 = ROOT / "p2040_s990_strict_same_scheme_subtraction_compatibility_residual_audit.py"
OUT2040 = ROOT / "generated" / "p2040_s990_strict_same_scheme_subtraction_compatibility_residual_audit.json"


class TestP2040SubtractionCompatibilityResidualAudit(unittest.TestCase):
    def test_p2040_computes_residual_bound_and_keeps_c3_open(self) -> None:
        subprocess.run([sys.executable, str(P2037)], check=True)
        subprocess.run([sys.executable, str(P2038)], check=True)
        subprocess.run([sys.executable, str(P2039)], check=True)
        subprocess.run([sys.executable, str(P2040)], check=True)

        data = json.loads(OUT2040.read_text(encoding="utf-8"))
        self.assertEqual(data["schema_version"], "p2040_s990_v1")
        self.assertEqual(data["status"], "OPEN_PARTIAL_PROGRESS_WITH_TRACE")
        self.assertEqual(
            data["result_kind"],
            "PASS_SAME_SCHEME_SUBTRACTION_COMPATIBILITY_RESIDUAL_AUDIT_WITH_BOUND__C3_STILL_OPEN",
        )

        ra = data["residual_audit"]
        self.assertGreaterEqual(ra["before_candidate_linf"], 0.0)
        self.assertGreaterEqual(ra["after_candidate_linf"], 0.0)
        self.assertGreater(ra["propagated_uncertainty_bound_linf"], 0.0)
        self.assertGreaterEqual(ra["after_candidate_with_uncertainty_upper_bound_linf"], ra["after_candidate_linf"])

        checks = data["gatekeeper_checks"]
        self.assertTrue(checks["preconditions_ready"])
        self.assertTrue(checks["base_candidate_nonzero"])
        self.assertTrue(checks["before_residual_finite"])
        self.assertTrue(checks["after_residual_finite"])
        self.assertTrue(checks["uncertainty_bound_finite"])
        self.assertFalse(checks["c3_theorem_proven"])
        self.assertTrue(checks["no_background_globalization_claimed"])
        self.assertTrue(checks["no_tensor_component_claimed"])
        self.assertTrue(checks["no_toe_closure_claimed"])

        c3 = data["c3_gate_update"]
        self.assertEqual(c3["C3_transport_theorem"], "OPEN")
        self.assertEqual(c3["C3_discharge_status"], "NOT_DISCHARGED")


if __name__ == "__main__":
    unittest.main()
