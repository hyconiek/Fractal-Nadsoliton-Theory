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
]
OUT = ROOT / "generated" / "p2045_s995_strict_same_scheme_seed_norm_joint_robustness_audit.json"


class TestP2045SeedNormJointRobustnessAudit(unittest.TestCase):
    def test_p2045_exports_joint_table_and_rank_stability(self) -> None:
        for script in SCRIPTS:
            subprocess.run([sys.executable, str(script)], check=True)

        data = json.loads(OUT.read_text(encoding="utf-8"))
        self.assertEqual(data["schema_version"], "p2045_s995_v1")
        self.assertEqual(data["status"], "OPEN_PARTIAL_PROGRESS_WITH_TRACE")
        self.assertEqual(
            data["result_kind"],
            "PASS_SEED_NORM_JOINT_ROBUSTNESS_AUDIT_WITH_RANK_STABILITY__C3_STILL_OPEN",
        )

        table = data["seed_norm_worst_case_table"]
        self.assertGreater(len(table), 0)

        rank = data["ranking_stability"]
        self.assertGreater(len(rank["pairwise"]), 0)
        self.assertGreaterEqual(rank["mean_spearman"], -1.0)
        self.assertLessEqual(rank["mean_spearman"], 1.0)
        self.assertGreaterEqual(rank["mean_kendall_tau"], -1.0)
        self.assertLessEqual(rank["mean_kendall_tau"], 1.0)

        checks = data["gatekeeper_checks"]
        self.assertTrue(checks["preconditions_ready"])
        self.assertTrue(checks["table_nonempty"])
        self.assertTrue(checks["pairwise_nonempty"])
        self.assertTrue(checks["mean_spearman_finite"])
        self.assertTrue(checks["mean_kendall_finite"])
        self.assertFalse(checks["c3_theorem_proven"])

        c3 = data["c3_gate_update"]
        self.assertEqual(c3["C3_transport_theorem"], "OPEN")
        self.assertEqual(c3["C3_discharge_status"], "NOT_DISCHARGED")


if __name__ == "__main__":
    unittest.main()
