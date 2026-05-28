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


class TestP2321S1271StrictDihedralNonlinearInvariantAxisCandidateAuditProbe(unittest.TestCase):
    def test_packet(self) -> None:
        subprocess.run(
            [sys.executable, str(ROOT / "p2321_s1271_strict_dihedral_nonlinear_invariant_axis_candidate_audit_probe.py")],
            check=True,
        )
        data = json.loads(
            (G / "p2321_s1271_strict_dihedral_nonlinear_invariant_axis_candidate_audit_probe.json").read_text(
                encoding="utf-8"
            )
        )
        self.assertEqual(data["schema_version"], "p2321_s1271_v1")
        self.assertEqual(data["packet_id"], "P2321")
        self.assertEqual(data["result_kind"], "STRICT_DIHEDRAL_NONLINEAR_INVARIANT_AXIS_CANDIDATE_AUDIT_NO_G1_G3_UPDATE")
        probe = data["strict_dihedral_nonlinear_invariant_axis_candidate_audit_probe"]
        self.assertGreaterEqual(probe["far_nonlinear_invariant_grep_audit"]["hit_count"], 5)
        rows = probe["pair_nonlinear_axis_rows"]
        self.assertEqual(len(rows), 5)
        expected_degrees = {"pair1": 12, "pair2": 6, "pair3": 4, "pair4": 3, "pair5": 12}
        for row in rows:
            self.assertEqual(row["lowest_d12_invariant_axis_degree"], expected_degrees[row["pair"]])
            self.assertLess(row["d12_sample_invariance_residual"], 1e-9)
            self.assertGreater(row["unit_circle_angular_spread"], 1.9)
            self.assertFalse(row["is_linear_or_quadratic"])
            self.assertFalse(row["supplies_signed_policy_orientation"])
            self.assertEqual(row["directed_extremal_ray_count"], 2 * row["rotation_order_q"])
            self.assertEqual(row["maxima_count"], row["rotation_order_q"])
            self.assertEqual(row["minima_count"], row["rotation_order_q"])
        cert = probe["nonlinear_axis_certificate"]
        self.assertEqual(cert["degree_sequence_by_pair"], expected_degrees)
        self.assertTrue(cert["all_pairs_have_nonlinear_d12_axis_harmonic"])
        self.assertTrue(cert["all_axis_harmonics_nonradial"])
        self.assertTrue(cert["all_axis_harmonics_outside_p2319_linear_quadratic_scope"])
        self.assertTrue(cert["all_axis_harmonics_fail_signed_policy_orientation"])
        self.assertLess(cert["max_d12_sample_invariance_residual"], 1e-9)
        update = probe["bridge_obligation_update"]
        self.assertAlmostEqual(update["required_lift_per_step"], 0.0068)
        self.assertTrue(update["nonlinear_axis_candidates_exported"])
        self.assertFalse(update["nonlinear_axis_candidates_fill_any_missing_p2318_field"])
        self.assertEqual(len(update["fields_still_unfilled_after_nonlinear_axis_candidates"]), 6)
        self.assertFalse(update["g1_g3_update_allowed"])
        theorem = probe["theorem_export"]
        self.assertEqual(sha256_json(theorem), probe["theorem_fingerprint_sha256"])
        self.assertIn("lowest-degree nonradial D12-invariant polynomial", theorem["formal_statement"])
        g = data["gatekeeper_checks"]
        self.assertTrue(g["far_nonlinear_invariant_grep_hits_found"])
        self.assertTrue(g["p2320_loaded"])
        self.assertTrue(g["p2318_loaded"])
        self.assertTrue(g["degree_sequence_matches_d12_pair_orders"])
        self.assertTrue(g["d12_invariance_residual_small"])
        self.assertTrue(g["all_axis_harmonics_nonradial"])
        self.assertTrue(g["outside_p2319_linear_quadratic_scope"])
        self.assertTrue(g["nonlinear_axis_candidates_not_promoted_to_bridge"])
        self.assertTrue(g["p2318_bridge_fields_still_unfilled"])
        self.assertTrue(g["strict_g1_g3_not_updated"])
        self.assertTrue(g["no_selector_closure_claimed"])
        self.assertTrue(g["no_qw2191_discharge_claimed"])
        self.assertTrue(g["no_legacy_kernel_role_transfer"])
        self.assertTrue(g["no_toe_closure_claimed"])


if __name__ == "__main__":
    unittest.main()
