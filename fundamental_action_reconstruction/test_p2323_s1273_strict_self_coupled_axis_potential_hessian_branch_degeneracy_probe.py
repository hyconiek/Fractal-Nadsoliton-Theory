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


class TestP2323S1273StrictSelfCoupledAxisPotentialHessianBranchDegeneracyProbe(unittest.TestCase):
    def test_packet(self) -> None:
        subprocess.run(
            [sys.executable, str(ROOT / "p2323_s1273_strict_self_coupled_axis_potential_hessian_branch_degeneracy_probe.py")],
            check=True,
        )
        data = json.loads(
            (G / "p2323_s1273_strict_self_coupled_axis_potential_hessian_branch_degeneracy_probe.json").read_text(
                encoding="utf-8"
            )
        )
        self.assertEqual(data["schema_version"], "p2323_s1273_v1")
        self.assertEqual(data["packet_id"], "P2323")
        self.assertEqual(data["result_kind"], "STRICT_SELF_COUPLED_AXIS_POTENTIAL_HESSIAN_BRANCH_DEGENERACY_AUDIT_NO_G1_G3_UPDATE")
        probe = data["strict_self_coupled_axis_potential_hessian_branch_degeneracy_probe"]
        self.assertGreaterEqual(probe["far_hessian_branch_degeneracy_grep_audit"]["hit_count"], 5)
        rows = probe["hessian_branch_rows"]
        self.assertEqual(len(rows), 5)
        expected_counts = {"pair1": 12, "pair2": 6, "pair3": 4, "pair4": 3, "pair5": 12}
        for row in rows:
            self.assertEqual(row["branch_count"], expected_counts[row["pair"]])
            self.assertTrue(row["all_branches_locally_stable"])
            self.assertTrue(row["all_branches_energy_degenerate"])
            self.assertLess(row["energy_spread_across_minima"], 1e-9)
            self.assertFalse(row["local_stability_supplies_branch_selector"])
            values = row["derivative_values"]
            self.assertGreater(values["radial_second_variation"], 0.0)
            self.assertGreater(values["angular_second_variation_arc_coordinate"], 0.0)
            self.assertLess(abs(values["radial_minus_arc_angular_second_variation"]), 1e-9)
            self.assertEqual(values["cross_second_variation_at_minimum"], 0.0)
            self.assertGreater(values["physical_hessian_determinant"], 0.0)
            self.assertLess(values["angular_second_variation_at_adjacent_maximum"], 0.0)
            for branch in row["branch_rows"]:
                self.assertTrue(branch["locally_stable_minimum"])
                self.assertEqual(len(branch["physical_hessian_eigenvalues"]), 2)
                self.assertTrue(all(eigenvalue > 0 for eigenvalue in branch["physical_hessian_eigenvalues"]))
        cert = probe["hessian_certificate"]
        self.assertEqual(cert["branch_count_by_pair"], expected_counts)
        self.assertTrue(cert["all_axis_branches_locally_stable"])
        self.assertTrue(cert["all_axis_branches_energy_degenerate"])
        self.assertLess(cert["max_energy_spread_across_minima"], 1e-9)
        self.assertGreater(cert["min_physical_hessian_eigenvalue"], 0.0)
        self.assertLess(cert["max_radial_arc_angular_hessian_difference"], 1e-9)
        self.assertFalse(cert["local_stability_breaks_branch_degeneracy"])
        update = probe["bridge_obligation_update"]
        self.assertAlmostEqual(update["required_lift_per_step"], 0.0068)
        self.assertTrue(update["hessian_stability_certificate_exported"])
        self.assertFalse(update["hessian_stability_fills_any_missing_p2318_field"])
        self.assertEqual(len(update["fields_still_unfilled_after_hessian_stability"]), 6)
        self.assertFalse(update["g1_g3_update_allowed"])
        theorem = probe["theorem_export"]
        self.assertEqual(sha256_json(theorem), probe["theorem_fingerprint_sha256"])
        self.assertIn("positive radial and physical angular second variations", theorem["formal_statement"])
        g = data["gatekeeper_checks"]
        self.assertTrue(g["far_hessian_branch_degeneracy_grep_hits_found"])
        self.assertTrue(g["p2322_loaded"])
        self.assertTrue(g["p2318_loaded"])
        self.assertTrue(g["all_axis_branches_locally_stable"])
        self.assertTrue(g["all_axis_branches_energy_degenerate"])
        self.assertTrue(g["physical_hessian_positive"])
        self.assertTrue(g["radial_and_arc_angular_hessians_match"])
        self.assertTrue(g["adjacent_maxima_detected_by_negative_angular_hessian"])
        self.assertTrue(g["hessian_stability_not_promoted_to_bridge"])
        self.assertTrue(g["p2318_bridge_fields_still_unfilled"])
        self.assertTrue(g["strict_g1_g3_not_updated"])
        self.assertTrue(g["no_selector_closure_claimed"])
        self.assertTrue(g["no_qw2191_discharge_claimed"])
        self.assertTrue(g["no_legacy_kernel_role_transfer"])
        self.assertTrue(g["no_toe_closure_claimed"])


if __name__ == "__main__":
    unittest.main()
