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
    ROOT / "p2058_s1008_strict_same_scheme_policy_regime_switch_audit.py",
    ROOT / "p2059_s1009_strict_same_scheme_robust_regime_regret_envelope_audit.py",
    ROOT / "p2060_s1010_strict_same_scheme_regime_aware_budget_reallocation_audit.py",
    ROOT / "p2061_s1011_strict_same_scheme_reallocation_stress_audit.py",
    ROOT / "p2062_s1012_strict_same_scheme_nonlinearity_stress_audit.py",
    ROOT / "p2063_s1013_strict_same_scheme_cross_model_robustness_certificate_audit.py",
    ROOT / "p2064_s1014_strict_same_scheme_certificate_perturbation_radius_audit.py",
    ROOT / "p2065_s1015_strict_same_scheme_certificate_transportability_gap_audit.py",
    ROOT / "p2067_s1017_strict_same_scheme_certificate_transportability_targeting_audit.py",
]
OUT = ROOT / "generated" / "p2067_s1017_strict_same_scheme_certificate_transportability_targeting_audit.json"


class TestP2067CertificateTransportabilityTargetingAudit(unittest.TestCase):
    def test_p2067_exports_targeting_and_keeps_c3_open(self) -> None:
        for script in SCRIPTS:
            subprocess.run([sys.executable, str(script)], check=True)

        data = json.loads(OUT.read_text(encoding="utf-8"))
        self.assertEqual(data["schema_version"], "p2067_s1017_v1")
        self.assertEqual(data["status"], "OPEN_PARTIAL_PROGRESS_WITH_TRACE")
        self.assertEqual(
            data["result_kind"],
            "PASS_CERTIFICATE_TRANSPORTABILITY_TARGETING_AUDIT_WITH_TRACE__C3_STILL_OPEN",
        )

        targeting = data["transportability_targeting"]
        self.assertGreaterEqual(targeting["transportability_gap"], 0.0)
        self.assertGreaterEqual(targeting["required_reduction_factor"], 0.0)
        self.assertLessEqual(targeting["required_reduction_factor"], 1.0)
        self.assertGreaterEqual(targeting["required_reduction_percent"], 0.0)

        checks = data["gatekeeper_checks"]
        self.assertTrue(checks["preconditions_ready"])
        self.assertTrue(checks["gap_nonnegative"])
        self.assertTrue(checks["required_reduction_factor_in_0_1"])
        self.assertTrue(checks["required_reduction_percent_nonnegative"])
        self.assertFalse(checks["c3_theorem_proven"])
        self.assertTrue(checks["no_background_globalization_claimed"])
        self.assertTrue(checks["no_tensor_component_claimed"])
        self.assertTrue(checks["no_toe_closure_claimed"])


if __name__ == "__main__":
    unittest.main()
