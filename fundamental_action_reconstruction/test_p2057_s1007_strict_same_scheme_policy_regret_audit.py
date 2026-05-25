from __future__ import annotations

import json
import subprocess
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent
SCRIPTS = [
    ROOT / "p2037_s987_strict_task1_same_scheme_finite_part_map_seed.py",
    ROOT / "p2038_s988_strict_same_scheme_finite_part_map_data_import_audit.py",
    ROOT / "p2039_s989_strict_same_scheme_finite_part_candidate_uncertainty_bound_probe.py",
    ROOT / "p2040_s990_strict_same_scheme_subtraction_compatibility_residual_audit.py",
    ROOT / "p2041_s991_strict_same_scheme_operator_level_subtraction_identity_witness.py",
    ROOT / "p2042_s992_strict_same_scheme_operator_level_residual_envelope_sweep.py",
    ROOT / "p2043_s993_strict_same_scheme_stratified_residual_envelope_audit.py",
    ROOT / "p2044_s994_strict_same_scheme_seed_sensitivity_meta_envelope_audit.py",
    ROOT / "p2045_s995_strict_same_scheme_seed_norm_joint_robustness_audit.py",
    ROOT / "p2046_s996_strict_same_scheme_adversarial_bucket_perturbation_audit.py",
    ROOT / "p2047_s997_strict_same_scheme_stability_certificate_extraction_audit.py",
    ROOT / "p2048_s998_strict_same_scheme_certificate_falsification_frontier_audit.py",
    ROOT / "p2049_s999_strict_same_scheme_frontier_confidence_audit.py",
    ROOT / "p2050_s1000_strict_same_scheme_censoring_resolution_audit.py",
    ROOT / "p2051_s1001_strict_same_scheme_budget_scaling_audit.py",
    ROOT / "p2052_s1002_strict_same_scheme_budget_normalized_frontier_efficiency_audit.py",
    ROOT / "p2053_s1003_strict_same_scheme_censoring_aware_budget_profile_recommendation_audit.py",
    ROOT / "p2054_s1004_strict_same_scheme_adaptive_budget_allocator_audit.py",
    ROOT / "p2055_s1005_strict_same_scheme_adaptive_vs_fixed_realized_frontier_gain_audit.py",
    ROOT / "p2056_s1006_strict_same_scheme_allocation_policy_calibration_audit.py",
    ROOT / "p2057_s1007_strict_same_scheme_policy_regret_audit.py",
]
OUT = ROOT / "generated" / "p2057_s1007_strict_same_scheme_policy_regret_audit.json"


class TestP2057PolicyRegretAudit(unittest.TestCase):
    def test_p2057_exports_regret_and_keeps_c3_open(self) -> None:
        for script in SCRIPTS:
            subprocess.run([sys.executable, str(script)], check=True)

        data = json.loads(OUT.read_text(encoding="utf-8"))
        self.assertEqual(data["schema_version"], "p2057_s1007_v1")
        self.assertEqual(data["status"], "OPEN_PARTIAL_PROGRESS_WITH_TRACE")
        self.assertEqual(data["result_kind"], "PASS_POLICY_REGRET_AUDIT_WITH_BOOTSTRAP_CI__C3_STILL_OPEN")

        self.assertIsNotNone(data["oracle_policy"])
        self.assertIsNotNone(data["selected_policy"])
        rs = data["regret_summary"]
        self.assertGreaterEqual(rs["regret_ci95"]["high"], rs["regret_ci95"]["low"])

        checks = data["gatekeeper_checks"]
        self.assertTrue(checks["preconditions_ready"])
        self.assertTrue(checks["policies_nonempty"])
        self.assertTrue(checks["oracle_present"])
        self.assertFalse(checks["c3_theorem_proven"])
        self.assertTrue(checks["no_background_globalization_claimed"])
        self.assertTrue(checks["no_tensor_component_claimed"])
        self.assertTrue(checks["no_toe_closure_claimed"])

        c3 = data["c3_gate_update"]
        self.assertEqual(c3["C3_transport_theorem"], "OPEN")
        self.assertEqual(c3["C3_discharge_status"], "NOT_DISCHARGED")


if __name__ == "__main__":
    unittest.main()
