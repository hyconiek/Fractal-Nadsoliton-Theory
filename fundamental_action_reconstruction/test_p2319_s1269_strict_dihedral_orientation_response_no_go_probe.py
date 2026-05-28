from __future__ import annotations

import hashlib
import json
import subprocess
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent
G = ROOT / "generated"


def sha256_json(payload: object) -> str:
    blob = json.dumps(payload, sort_keys=True, separators=(",", ":"), ensure_ascii=False)
    return hashlib.sha256(blob.encode("utf-8")).hexdigest()


class TestP2319S1269StrictDihedralOrientationResponseNoGoProbe(unittest.TestCase):
    def test_packet(self) -> None:
        subprocess.run(
            [sys.executable, str(ROOT / "p2319_s1269_strict_dihedral_orientation_response_no_go_probe.py")],
            check=True,
        )
        data = json.loads(
            (G / "p2319_s1269_strict_dihedral_orientation_response_no_go_probe.json").read_text(
                encoding="utf-8"
            )
        )
        self.assertEqual(data["schema_version"], "p2319_s1269_v1")
        self.assertEqual(data["packet_id"], "P2319")
        self.assertEqual(data["result_kind"], "STRICT_DIHEDRAL_ORIENTATION_RESPONSE_NO_GO_NO_G1_G3_UPDATE")
        probe = data["strict_dihedral_orientation_response_no_go_probe"]
        self.assertGreaterEqual(probe["far_symmetry_response_grep_audit"]["hit_count"], 5)
        self.assertEqual(probe["d12_group_order"], 24)
        rows = probe["pair_orientation_reynolds_rows"]
        self.assertEqual(len(rows), 5)
        for row in rows:
            self.assertLess(row["kernel_scalar_identity_residual"], 1e-10)
            self.assertLess(row["reynolds_linear_cos_norm"], 1e-10)
            self.assertLess(row["reynolds_linear_sin_norm"], 1e-10)
            self.assertLess(row["reynolds_traceless_anisotropy_inf_norm"], 1e-10)
            self.assertLess(row["reynolds_shear_inf_norm"], 1e-10)
            self.assertLess(row["reynolds_handed_skew_inf_norm"], 1e-10)
            self.assertLess(row["projector_reynolds_preservation_inf_residual"], 1e-10)
            self.assertFalse(row["orientation_odd_linear_survives_d12_average"])
            self.assertFalse(row["orientation_anisotropy_survives_d12_average"])
            self.assertFalse(row["handed_pseudoscalar_survives_d12_average"])
            self.assertTrue(row["only_unsigned_pair_norm_survives"])
        cert = probe["no_go_certificate"]
        self.assertLess(cert["kernel_d12_invariance_residual"], 1e-10)
        self.assertLess(cert["alpha_scaled_kernel_d12_invariance_residual"], 1e-10)
        self.assertTrue(cert["all_orientation_odd_linear_reynolds_projections_zero"])
        self.assertTrue(cert["all_traceless_quadratic_orientation_projections_zero"])
        self.assertTrue(cert["all_handed_skew_projections_zero_under_reflection"])
        self.assertTrue(cert["all_unsigned_pair_norm_projectors_survive"])
        update = probe["bridge_obligation_update"]
        self.assertAlmostEqual(update["required_lift_per_step"], 0.0068)
        self.assertEqual(len(update["fields_still_unfilled_by_d12_kernel_scalar_class"]), 6)
        self.assertFalse(update["d12_class_supplies_signed_scalar_lift_per_step"])
        self.assertFalse(update["d12_class_supplies_margin_response_functional"])
        self.assertFalse(update["g1_g3_update_allowed"])
        theorem = probe["theorem_export"]
        self.assertEqual(sha256_json(theorem), probe["theorem_fingerprint_sha256"])
        self.assertIn("no signed Fourier-pair orientation response survives", theorem["formal_statement"])
        g = data["gatekeeper_checks"]
        self.assertTrue(g["far_symmetry_response_grep_hits_found"])
        self.assertTrue(g["p2317_loaded"])
        self.assertTrue(g["p2318_loaded"])
        self.assertTrue(g["p2308_current_class_nonexistence_loaded"])
        self.assertTrue(g["d12_group_order_is_24"])
        self.assertTrue(g["kernel_d12_invariance_verified"])
        self.assertTrue(g["alpha_scaled_kernel_d12_invariance_verified"])
        self.assertTrue(g["orientation_odd_linear_responses_project_to_zero"])
        self.assertTrue(g["traceless_quadratic_orientation_responses_project_to_zero"])
        self.assertTrue(g["handed_reflection_odd_responses_project_to_zero"])
        self.assertTrue(g["only_unsigned_pair_norms_survive"])
        self.assertTrue(g["p2318_bridge_fields_still_unfilled"])
        self.assertTrue(g["strict_g1_g3_not_updated"])
        self.assertTrue(g["no_selector_closure_claimed"])
        self.assertTrue(g["no_qw2191_discharge_claimed"])
        self.assertTrue(g["no_legacy_kernel_role_transfer"])
        self.assertTrue(g["no_toe_closure_claimed"])


if __name__ == "__main__":
    unittest.main()
