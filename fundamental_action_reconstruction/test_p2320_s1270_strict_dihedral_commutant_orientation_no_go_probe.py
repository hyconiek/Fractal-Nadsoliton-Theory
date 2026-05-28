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


class TestP2320S1270StrictDihedralCommutantOrientationNoGoProbe(unittest.TestCase):
    def test_packet(self) -> None:
        subprocess.run(
            [sys.executable, str(ROOT / "p2320_s1270_strict_dihedral_commutant_orientation_no_go_probe.py")],
            check=True,
        )
        data = json.loads(
            (G / "p2320_s1270_strict_dihedral_commutant_orientation_no_go_probe.json").read_text(
                encoding="utf-8"
            )
        )
        self.assertEqual(data["schema_version"], "p2320_s1270_v1")
        self.assertEqual(data["packet_id"], "P2320")
        self.assertEqual(data["result_kind"], "STRICT_DIHEDRAL_COMMUTANT_ORIENTATION_NO_GO_NO_G1_G3_UPDATE")
        probe = data["strict_dihedral_commutant_orientation_no_go_probe"]
        self.assertGreaterEqual(probe["far_commutant_grep_audit"]["hit_count"], 5)
        cert = probe["commutant_certificate"]
        self.assertEqual(cert["d12_group_order"], 24)
        self.assertEqual(cert["c12_group_order"], 12)
        self.assertEqual(cert["d12_distance_basis_count"], 7)
        self.assertEqual(cert["d12_distance_basis_gram_rank"], 7)
        self.assertEqual(cert["c12_shift_basis_count"], 12)
        self.assertEqual(cert["c12_shift_basis_gram_rank"], 12)
        self.assertTrue(cert["all_distance_basis_d12_invariant"])
        self.assertTrue(cert["all_shift_basis_c12_invariant"])
        self.assertTrue(cert["some_shift_basis_reflection_noninvariant"])
        self.assertTrue(cert["all_non_unsigned_orientation_candidates_project_to_zero_in_d12_commutant"])
        self.assertTrue(cert["all_unsigned_pair_projectors_survive_d12_commutant"])
        self.assertEqual(cert["c12_handed_skew_survivor_count"], 5)
        self.assertTrue(cert["c12_handed_skew_survivors_are_reflection_odd_not_d12_admissible"])
        self.assertLess(cert["max_d12_non_unsigned_projection_norm"], 1e-10)
        self.assertGreater(cert["min_d12_unsigned_projection_norm"], 0.9)
        self.assertEqual(len(probe["d12_basis_rows"]), 7)
        self.assertEqual(len(probe["c12_contrast_basis_rows"]), 12)
        self.assertEqual(len(probe["pair_operator_projection_rows"]), 5)
        for pair in probe["pair_operator_projection_rows"]:
            candidates = {row["candidate_id"]: row for row in pair["candidate_rows"]}
            self.assertTrue(candidates["unsigned_pair_projector"]["d12_survives_commutant_projection"])
            self.assertFalse(candidates["traceless_anisotropy_cc_minus_ss"]["d12_survives_commutant_projection"])
            self.assertFalse(candidates["symmetric_shear_cs_plus_sc"]["d12_survives_commutant_projection"])
            self.assertFalse(candidates["antisymmetric_handed_skew_cs_minus_sc"]["d12_survives_commutant_projection"])
            self.assertTrue(candidates["antisymmetric_handed_skew_cs_minus_sc"]["c12_survives_commutant_projection"])
        update = probe["bridge_obligation_update"]
        self.assertAlmostEqual(update["required_lift_per_step"], 0.0068)
        self.assertFalse(update["commutant_class_fills_any_missing_p2318_field"])
        self.assertEqual(len(update["fields_still_unfilled_by_full_d12_commutant_class"]), 6)
        self.assertFalse(update["g1_g3_update_allowed"])
        theorem = probe["theorem_export"]
        self.assertEqual(sha256_json(theorem), probe["theorem_fingerprint_sha256"])
        self.assertIn("full D12-invariant linear-operator commutant", theorem["formal_statement"])
        g = data["gatekeeper_checks"]
        self.assertTrue(g["far_commutant_grep_hits_found"])
        self.assertTrue(g["p2319_loaded"])
        self.assertTrue(g["p2318_loaded"])
        self.assertTrue(g["d12_commutant_dimension_is_7"])
        self.assertTrue(g["c12_commutant_dimension_is_12"])
        self.assertTrue(g["non_unsigned_orientation_operators_annihilated_by_d12_commutant_projection"])
        self.assertTrue(g["unsigned_pair_projectors_survive_d12_commutant_projection"])
        self.assertTrue(g["c12_handed_contrast_detected_but_not_promoted"])
        self.assertTrue(g["p2318_bridge_fields_still_unfilled"])
        self.assertTrue(g["strict_g1_g3_not_updated"])
        self.assertTrue(g["no_selector_closure_claimed"])
        self.assertTrue(g["no_qw2191_discharge_claimed"])
        self.assertTrue(g["no_legacy_kernel_role_transfer"])
        self.assertTrue(g["no_toe_closure_claimed"])


if __name__ == "__main__":
    unittest.main()
