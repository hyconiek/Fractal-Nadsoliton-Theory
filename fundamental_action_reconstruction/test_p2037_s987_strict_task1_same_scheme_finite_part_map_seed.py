from __future__ import annotations

import json
import subprocess
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent
P2034 = ROOT / "p2034_s984_strict_task1_quotient_only_renormalization_theorem.py"
P2035 = ROOT / "p2035_s985_strict_task1_quotient_background_transport_obstruction_theorem.py"
P2036 = ROOT / "p2036_s986_strict_task1_quotient_background_transport_candidate_contract.py"
P2037 = ROOT / "p2037_s987_strict_task1_same_scheme_finite_part_map_seed.py"
OUT2037 = ROOT / "generated" / "p2037_s987_strict_task1_same_scheme_finite_part_map_seed.json"


class TestP2037SameSchemeFinitePartMapSeed(unittest.TestCase):
    def test_p2037_exports_first_explicit_same_scheme_seed_without_false_discharge(self) -> None:
        # Keep this test local to P2037: dependencies are read-only inputs and
        # should not be regenerated here to avoid unrelated artifact churn.
        for dependency in (P2034, P2035, P2036):
            self.assertTrue(dependency.exists(), f"Missing dependency script: {dependency}")
        subprocess.run([sys.executable, str(P2037)], check=True)

        data = json.loads(OUT2037.read_text(encoding="utf-8"))

        self.assertEqual(data["schema_version"], "p2037_s987_v1")
        self.assertEqual(data["status"], "OPEN_PARTIAL_PROGRESS_WITH_TRACE")
        self.assertEqual(
            data["result_kind"],
            "PASS_EXPLICIT_SAME_SCHEME_FINITE_PART_MAP_SEED_WITH_TRACE__C2_C3_STILL_OPEN",
        )

        seed = data["same_scheme_finite_part_map_seed"]
        self.assertEqual(seed["map_id"], "M_same_scheme_seed_v1")
        self.assertEqual(seed["basis"], ["R2_bar", "Ric2_bar", "Riem2_bar"])
        self.assertEqual(seed["sigma2_anchor"], 0.0)
        self.assertEqual(seed["delta_finite_part_matrix"], [[0.0, 0.0, 0.0], [0.0, 0.0, 0.0], [0.0, 0.0, 0.0]])
        self.assertEqual(seed["transport_matrix"], [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]])
        self.assertEqual(seed["residual_vs_identity_linf"], 0.0)

        c2c3 = data["c2_c3_gate_update"]
        self.assertEqual(c2c3["C2_basis_preserving_map_seed"], "EXPORTED")
        self.assertEqual(c2c3["C3_finite_part_scheme_transport_theorem"], "OPEN")
        self.assertEqual(c2c3["C2_C3_discharge_status"], "NOT_DISCHARGED")

        checks = data["gatekeeper_checks"]
        self.assertTrue(checks["map_seed_exported"])
        self.assertTrue(checks["residual_vs_identity_zero"])
        self.assertTrue(checks["c2_basis_preserving_seed_exported"])
        self.assertFalse(checks["c3_theorem_proven"])
        self.assertTrue(checks["no_background_globalization_claimed"])
        self.assertTrue(checks["no_tensor_component_claimed"])
        self.assertTrue(checks["no_toe_closure_claimed"])


if __name__ == "__main__":
    unittest.main()
