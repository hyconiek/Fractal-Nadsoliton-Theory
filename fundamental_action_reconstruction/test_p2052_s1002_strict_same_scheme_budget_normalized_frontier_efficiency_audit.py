from __future__ import annotations

import json
import math
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
]
OUT = ROOT / "generated" / "p2052_s1002_strict_same_scheme_budget_normalized_frontier_efficiency_audit.json"


class TestP2052BudgetNormalizedFrontierEfficiencyAudit(unittest.TestCase):
    def test_p2052_exports_efficiency_rows_and_keeps_c3_open(self) -> None:
        for script in SCRIPTS:
            subprocess.run([sys.executable, str(script)], check=True)

        data = json.loads(OUT.read_text(encoding="utf-8"))
        self.assertEqual(data["schema_version"], "p2052_s1002_v1")
        self.assertEqual(data["status"], "OPEN_PARTIAL_PROGRESS_WITH_TRACE")
        self.assertEqual(
            data["result_kind"],
            "PASS_BUDGET_NORMALIZED_FRONTIER_EFFICIENCY_AUDIT_WITH_TRACE__C3_STILL_OPEN",
        )

        rows = data["profile_efficiency_rows"]
        self.assertGreater(len(rows), 0)

        for row in rows:
            self.assertGreaterEqual(row["scan_steps_used"], 0)
            self.assertGreaterEqual(row["final_sign_pattern_budget"], 0)
            self.assertGreaterEqual(row["vectors_per_bucket"], 0)
            self.assertGreaterEqual(row["compute_units_proxy"], 1)
            self.assertIn(row["detected_frontier_flag"], (0, 1))
            self.assertGreaterEqual(row["frontier_detection_per_compute_unit"], 0.0)

        summary = data["efficiency_summary"]
        self.assertGreaterEqual(summary["detected_frontier_count"], 0)
        self.assertLessEqual(summary["detected_frontier_count"], summary["profiles_total"])
        self.assertTrue(math.isfinite(summary["mean_frontier_detection_per_compute_unit"]))

        checks = data["gatekeeper_checks"]
        self.assertTrue(checks["preconditions_ready"])
        self.assertTrue(checks["rows_nonempty"])
        self.assertTrue(checks["detected_count_consistent"])
        self.assertTrue(checks["mean_eff_finite"])
        self.assertFalse(checks["c3_theorem_proven"])
        self.assertTrue(checks["no_background_globalization_claimed"])
        self.assertTrue(checks["no_tensor_component_claimed"])
        self.assertTrue(checks["no_toe_closure_claimed"])

        c3 = data["c3_gate_update"]
        self.assertEqual(c3["C3_transport_theorem"], "OPEN")
        self.assertEqual(c3["C3_discharge_status"], "NOT_DISCHARGED")


if __name__ == "__main__":
    unittest.main()
