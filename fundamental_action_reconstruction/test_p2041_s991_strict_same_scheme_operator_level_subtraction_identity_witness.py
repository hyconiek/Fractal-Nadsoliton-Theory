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
P2041 = ROOT / "p2041_s991_strict_same_scheme_operator_level_subtraction_identity_witness.py"
OUT2041 = ROOT / "generated" / "p2041_s991_strict_same_scheme_operator_level_subtraction_identity_witness.json"


class TestP2041OperatorLevelSubtractionIdentityWitness(unittest.TestCase):
    def test_p2041_exports_symbolic_numeric_witness_and_keeps_c3_open(self) -> None:
        subprocess.run([sys.executable, str(P2037)], check=True)
        subprocess.run([sys.executable, str(P2038)], check=True)
        subprocess.run([sys.executable, str(P2039)], check=True)
        subprocess.run([sys.executable, str(P2040)], check=True)
        subprocess.run([sys.executable, str(P2041)], check=True)

        data = json.loads(OUT2041.read_text(encoding="utf-8"))
        self.assertEqual(data["schema_version"], "p2041_s991_v1")
        self.assertEqual(data["status"], "OPEN_PARTIAL_PROGRESS_WITH_TRACE")
        self.assertEqual(
            data["result_kind"],
            "PASS_OPERATOR_LEVEL_SAME_SCHEME_SUBTRACTION_IDENTITY_WITNESS_WITH_BOUND__C3_STILL_OPEN",
        )

        sw = data["symbolic_witness"]
        self.assertIn("-delta*x", sw["identity_form"])
        self.assertEqual(len(sw["residual_symbolic_components"]), 3)

        nw = data["numeric_witness"]
        self.assertEqual(nw["test_vectors_count"], 4)
        self.assertGreaterEqual(nw["worst_case_residual_linf"], 0.0)
        self.assertGreaterEqual(nw["worst_case_residual_with_uncertainty_bound_linf"], nw["worst_case_residual_linf"])

        checks = data["gatekeeper_checks"]
        self.assertTrue(checks["preconditions_ready"])
        self.assertTrue(checks["test_vectors_count_ge_4"])
        self.assertTrue(checks["worst_case_finite"])
        self.assertTrue(checks["symbolic_identity_exported"])
        self.assertFalse(checks["c3_theorem_proven"])
        self.assertTrue(checks["no_background_globalization_claimed"])
        self.assertTrue(checks["no_tensor_component_claimed"])
        self.assertTrue(checks["no_toe_closure_claimed"])

        c3 = data["c3_gate_update"]
        self.assertEqual(c3["C3_transport_theorem"], "OPEN")
        self.assertEqual(c3["C3_discharge_status"], "NOT_DISCHARGED")


if __name__ == "__main__":
    unittest.main()
