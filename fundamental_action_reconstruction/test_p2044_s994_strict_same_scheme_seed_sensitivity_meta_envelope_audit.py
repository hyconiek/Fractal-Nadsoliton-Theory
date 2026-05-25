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
P2042 = ROOT / "p2042_s992_strict_same_scheme_operator_level_residual_envelope_sweep.py"
P2043 = ROOT / "p2043_s993_strict_same_scheme_stratified_residual_envelope_audit.py"
P2044 = ROOT / "p2044_s994_strict_same_scheme_seed_sensitivity_meta_envelope_audit.py"
OUT2044 = ROOT / "generated" / "p2044_s994_strict_same_scheme_seed_sensitivity_meta_envelope_audit.json"


class TestP2044SeedSensitivityMetaEnvelopeAudit(unittest.TestCase):
    def test_p2044_exports_meta_envelope_and_keeps_c3_open(self) -> None:
        for s in (P2037, P2038, P2039, P2040, P2041, P2042, P2043, P2044):
            subprocess.run([sys.executable, str(s)], check=True)

        data = json.loads(OUT2044.read_text(encoding="utf-8"))
        self.assertEqual(data["schema_version"], "p2044_s994_v1")
        self.assertEqual(data["status"], "OPEN_PARTIAL_PROGRESS_WITH_TRACE")
        self.assertEqual(
            data["result_kind"],
            "PASS_SEED_SENSITIVITY_META_ENVELOPE_AUDIT_WITH_INTERSEED_CI__C3_STILL_OPEN",
        )

        spec = data["seed_sensitivity_spec"]
        self.assertGreaterEqual(len(spec["seeds"]), 5)

        meta = data["meta_envelope"]
        self.assertGreater(meta["meta_sup_residual_with_buffer_linf"], 0.0)
        self.assertLessEqual(meta["interseed_ci_low"], meta["interseed_ci_high"])

        checks = data["gatekeeper_checks"]
        self.assertTrue(checks["preconditions_ready"])
        self.assertTrue(checks["seed_count_ge_5"])
        self.assertTrue(checks["meta_sup_finite"])
        self.assertTrue(checks["interseed_ci_finite"])
        self.assertTrue(checks["nonzero_delta_present"])
        self.assertFalse(checks["c3_theorem_proven"])
        self.assertTrue(checks["no_background_globalization_claimed"])
        self.assertTrue(checks["no_tensor_component_claimed"])
        self.assertTrue(checks["no_toe_closure_claimed"])

        c3 = data["c3_gate_update"]
        self.assertEqual(c3["C3_transport_theorem"], "OPEN")
        self.assertEqual(c3["C3_discharge_status"], "NOT_DISCHARGED")


if __name__ == "__main__":
    unittest.main()
