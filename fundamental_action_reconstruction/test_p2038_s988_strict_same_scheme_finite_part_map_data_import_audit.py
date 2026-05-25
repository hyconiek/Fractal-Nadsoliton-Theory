from __future__ import annotations

import json
import subprocess
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent
P2037 = ROOT / "p2037_s987_strict_task1_same_scheme_finite_part_map_seed.py"
P2038 = ROOT / "p2038_s988_strict_same_scheme_finite_part_map_data_import_audit.py"
OUT2038 = ROOT / "generated" / "p2038_s988_strict_same_scheme_finite_part_map_data_import_audit.json"


class TestP2038SameSchemeFinitePartImportAudit(unittest.TestCase):
    def test_p2038_imports_first_nonzero_candidate_without_false_c3_discharge(self) -> None:
        subprocess.run([sys.executable, str(P2037)], check=True)
        subprocess.run([sys.executable, str(P2038)], check=True)

        data = json.loads(OUT2038.read_text(encoding="utf-8"))

        self.assertEqual(data["schema_version"], "p2038_s988_v1")
        self.assertEqual(data["status"], "OPEN_PARTIAL_PROGRESS_WITH_TRACE")
        self.assertEqual(
            data["result_kind"],
            "PASS_FIRST_NONZERO_SAME_SCHEME_FINITE_PART_CANDIDATE_IMPORTED_WITH_TRACE__C3_STILL_OPEN",
        )

        imported = data["candidate_data_import"]
        self.assertEqual(imported["basis"], ["R2_bar", "Ric2_bar", "Riem2_bar"])
        self.assertEqual(imported["controlled_background_pair"], "B1_scalar_to_FRW_controlled_pair_v1")
        self.assertGreater(imported["candidate_norm_linf"], 0.0)

        c2c3 = data["c2_c3_gate_update"]
        self.assertEqual(c2c3["C2_basis_preserving_seed"], "KEPT_FROM_P2037")
        self.assertEqual(c2c3["C3_nonzero_candidate_data"], "IMPORTED_FOR_ONE_CONTROLLED_PAIR")
        self.assertEqual(c2c3["C3_finite_part_scheme_transport_theorem"], "OPEN")
        self.assertEqual(c2c3["C3_discharge_status"], "NOT_DISCHARGED")

        checks = data["gatekeeper_checks"]
        self.assertTrue(checks["seed_ready"])
        self.assertTrue(checks["nonzero_candidate_present"])
        self.assertTrue(checks["candidate_norm_positive"])
        self.assertFalse(checks["c3_theorem_proven"])
        self.assertTrue(checks["no_background_globalization_claimed"])
        self.assertTrue(checks["no_tensor_component_claimed"])
        self.assertTrue(checks["no_toe_closure_claimed"])


if __name__ == "__main__":
    unittest.main()
