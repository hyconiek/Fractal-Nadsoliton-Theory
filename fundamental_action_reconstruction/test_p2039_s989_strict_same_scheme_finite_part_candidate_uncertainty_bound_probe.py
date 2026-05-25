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
OUT2039 = ROOT / "generated" / "p2039_s989_strict_same_scheme_finite_part_candidate_uncertainty_bound_probe.json"


class TestP2039UncertaintyBoundProbe(unittest.TestCase):
    def test_p2039_computes_first_uncertainty_bound_and_keeps_c3_open(self) -> None:
        subprocess.run([sys.executable, str(P2037)], check=True)
        subprocess.run([sys.executable, str(P2038)], check=True)
        subprocess.run([sys.executable, str(P2039)], check=True)

        data = json.loads(OUT2039.read_text(encoding="utf-8"))
        self.assertEqual(data["schema_version"], "p2039_s989_v1")
        self.assertEqual(data["status"], "OPEN_PARTIAL_PROGRESS_WITH_TRACE")
        self.assertEqual(
            data["result_kind"],
            "PASS_FIRST_COMPUTED_UNCERTAINTY_BOUND_FOR_NONZERO_FINITE_PART_CANDIDATE_WITH_TRACE__C3_STILL_OPEN",
        )

        ub = data["uncertainty_bound"]
        self.assertEqual(ub["bound_type"], "small_grid_linf_uncertainty_bound")
        self.assertEqual(ub["grid_size"], 81)
        self.assertGreater(ub["absolute_linf_bound"], 0.0)
        self.assertGreaterEqual(ub["relative_linf_bound_vs_base"], 0.0)

        checks = data["gatekeeper_checks"]
        self.assertTrue(checks["preconditions_ready"])
        self.assertTrue(checks["grid_size_81"])
        self.assertTrue(checks["base_candidate_nonzero"])
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
