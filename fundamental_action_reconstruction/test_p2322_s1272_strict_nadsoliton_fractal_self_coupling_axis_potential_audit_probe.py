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


class TestP2322S1272StrictNadsolitonFractalSelfCouplingAxisPotentialAuditProbe(unittest.TestCase):
    def test_packet(self) -> None:
        subprocess.run(
            [sys.executable, str(ROOT / "p2322_s1272_strict_nadsoliton_fractal_self_coupling_axis_potential_audit_probe.py")],
            check=True,
        )
        data = json.loads(
            (G / "p2322_s1272_strict_nadsoliton_fractal_self_coupling_axis_potential_audit_probe.json").read_text(
                encoding="utf-8"
            )
        )
        self.assertEqual(data["schema_version"], "p2322_s1272_v1")
        self.assertEqual(data["packet_id"], "P2322")
        self.assertEqual(data["result_kind"], "STRICT_NADSOLITON_FRACTAL_SELF_COUPLING_AXIS_POTENTIAL_AUDIT_NO_G1_G3_UPDATE")
        probe = data["strict_nadsoliton_fractal_self_coupling_axis_potential_audit_probe"]
        self.assertGreaterEqual(probe["far_nadsoliton_self_coupling_grep_audit"]["hit_count"], 5)
        ontology = probe["ontology_audit"]
        self.assertTrue(ontology["nadsoliton_as_sole_primordial_information_respected"])
        self.assertTrue(ontology["no_separate_informational_layer_introduced"])
        self.assertEqual(ontology["preferred_internal_order"], "nadsoliton -> light -> matter -> emergent observer")
        rows = probe["self_coupling_axis_potential_rows"]
        self.assertEqual(len(rows), 5)
        expected = {"pair1": 12, "pair2": 6, "pair3": 4, "pair4": 3, "pair5": 12}
        for row in rows:
            self.assertEqual(row["degree_q"], expected[row["pair"]])
            self.assertTrue(row["real_part_matches_p2321_coefficients"])
            self.assertGreaterEqual(row["self_coupling_multiplication_depth"], 2)
            potential = row["fractal_self_potential_candidate"]
            self.assertTrue(potential["uses_only_same_pair_amplitude_self_couplings"])
            self.assertTrue(potential["does_not_introduce_sub_nadsoliton_information_layer"])
            self.assertTrue(potential["bounded_below_for_lambda_positive"])
            self.assertEqual(potential["minima_count"], row["degree_q"])
            self.assertEqual(potential["maxima_count"], row["degree_q"])
            self.assertAlmostEqual(potential["critical_radius_power_r_to_q"], 0.5)
            self.assertAlmostEqual(potential["minimum_value"], -0.25)
            self.assertLess(potential["max_abs_radial_gradient_at_minima"], 1e-9)
            self.assertLess(potential["max_abs_angular_gradient_at_minima"], 1e-9)
        cert = probe["self_coupling_certificate"]
        self.assertEqual(cert["target_degrees"], [3, 4, 6, 12])
        self.assertTrue(cert["all_p2321_polynomials_reproduced_by_iterated_self_coupling"])
        self.assertEqual(cert["max_self_coupling_multiplication_depth"], 4)
        self.assertTrue(cert["all_potentials_bounded_below"])
        self.assertTrue(cert["all_axis_minima_counts_match_degree"])
        update = probe["bridge_obligation_update"]
        self.assertAlmostEqual(update["required_lift_per_step"], 0.0068)
        self.assertTrue(update["self_coupling_axis_potentials_exported"])
        self.assertFalse(update["self_coupling_axis_potentials_fill_any_missing_p2318_field"])
        self.assertEqual(len(update["fields_still_unfilled_after_self_coupling_potential_candidates"]), 6)
        self.assertFalse(update["g1_g3_update_allowed"])
        theorem = probe["theorem_export"]
        self.assertEqual(sha256_json(theorem), probe["theorem_fingerprint_sha256"])
        self.assertIn("iterated self-couplings of the same pair-plane nadsoliton amplitude", theorem["formal_statement"])
        g = data["gatekeeper_checks"]
        self.assertTrue(g["far_nadsoliton_self_coupling_grep_hits_found"])
        self.assertTrue(g["p2321_loaded"])
        self.assertTrue(g["p2318_loaded"])
        self.assertTrue(g["nadsoliton_as_sole_primordial_information_respected"])
        self.assertTrue(g["no_separate_informational_layer_introduced"])
        self.assertTrue(g["all_p2321_polynomials_reproduced_by_self_coupling"])
        self.assertTrue(g["all_potentials_bounded_below"])
        self.assertTrue(g["potential_gradient_residuals_small"])
        self.assertTrue(g["axis_minima_counts_match_degree"])
        self.assertTrue(g["self_coupling_candidates_not_promoted_to_bridge"])
        self.assertTrue(g["p2318_bridge_fields_still_unfilled"])
        self.assertTrue(g["strict_g1_g3_not_updated"])
        self.assertTrue(g["no_selector_closure_claimed"])
        self.assertTrue(g["no_qw2191_discharge_claimed"])
        self.assertTrue(g["no_legacy_kernel_role_transfer"])
        self.assertTrue(g["no_toe_closure_claimed"])


if __name__ == "__main__":
    unittest.main()
