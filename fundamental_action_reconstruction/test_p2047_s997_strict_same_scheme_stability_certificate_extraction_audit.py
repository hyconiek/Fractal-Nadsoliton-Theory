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
]
OUT = ROOT / "generated" / "p2047_s997_strict_same_scheme_stability_certificate_extraction_audit.json"


class TestP2047StabilityCertificateExtractionAudit(unittest.TestCase):
    def test_p2047_exports_local_certificate_and_keeps_c3_open(self) -> None:
        for script in SCRIPTS:
            subprocess.run([sys.executable, str(script)], check=True)

        data = json.loads(OUT.read_text(encoding="utf-8"))
        self.assertEqual(data["schema_version"], "p2047_s997_v1")
        self.assertEqual(data["status"], "OPEN_PARTIAL_PROGRESS_WITH_TRACE")
        self.assertEqual(
            data["result_kind"],
            "PASS_LOCAL_STABILITY_CERTIFICATE_EXTRACTED_WITH_TRACE__C3_STILL_OPEN",
        )

        cert = data["certificate_extraction"]
        self.assertGreaterEqual(cert["epsilon_star_max_stable"], 0.0)
        self.assertEqual(len(cert["epsilon_grid"]), len(cert["stable_flags"]))
        self.assertGreater(cert["adversarial_sign_patterns_checked"], 0)

        checks = data["gatekeeper_checks"]
        self.assertTrue(checks["preconditions_ready"])
        self.assertTrue(checks["epsilon_star_finite"])
        self.assertTrue(checks["epsilon_star_nonnegative"])
        self.assertTrue(checks["grid_nonempty"])
        self.assertTrue(checks["flags_match_grid"])
        self.assertFalse(checks["c3_theorem_proven"])
        self.assertTrue(checks["no_background_globalization_claimed"])
        self.assertTrue(checks["no_tensor_component_claimed"])
        self.assertTrue(checks["no_toe_closure_claimed"])

        c3 = data["c3_gate_update"]
        self.assertEqual(c3["C3_transport_theorem"], "OPEN")
        self.assertEqual(c3["C3_discharge_status"], "NOT_DISCHARGED")


if __name__ == "__main__":
    unittest.main()
