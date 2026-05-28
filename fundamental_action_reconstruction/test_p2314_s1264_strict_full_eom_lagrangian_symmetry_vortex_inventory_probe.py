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
    blob = json.dumps(payload, sort_keys=True, separators=(",", ":"))
    return hashlib.sha256(blob.encode("utf-8")).hexdigest()


class TestP2314S1264StrictFullEomLagrangianSymmetryVortexInventoryProbe(unittest.TestCase):
    def test_packet(self) -> None:
        subprocess.run(
            [sys.executable, str(ROOT / "p2314_s1264_strict_full_eom_lagrangian_symmetry_vortex_inventory_probe.py")],
            check=True,
        )
        data = json.loads(
            (G / "p2314_s1264_strict_full_eom_lagrangian_symmetry_vortex_inventory_probe.json").read_text(encoding="utf-8")
        )
        self.assertEqual(data["schema_version"], "p2314_s1264_v1")
        self.assertEqual(data["packet_id"], "P2314")
        probe = data["strict_full_eom_lagrangian_symmetry_vortex_inventory_probe"]
        eom = probe["eom_lagrangian_inventory"]
        self.assertTrue(eom["schematic_candidate_lagrangian_present"])
        self.assertTrue(eom["candidate_lagrangian_eom_layer_present_in_release7"])
        self.assertTrue(eom["partial_adm_bianchi_lapse_chain_present"])
        self.assertTrue(eom["partial_spatial_eom_provider_matrix_present"])
        self.assertTrue(eom["full_tensor_component_profile_table_missing"])
        self.assertTrue(eom["curved_metric_projection_rule_missing"])
        self.assertFalse(eom["full_strict_eom_exported"])
        self.assertFalse(eom["full_theorem_grade_lagrangian_for_task3_exported"])
        alpha = probe["alpha_4ln2_symmetry_audit"]
        self.assertTrue(alpha["alpha_geo_strict_source_loaded"])
        self.assertTrue(alpha["computed_entropy_equals_4_ln2"])
        self.assertTrue(alpha["permutation_invariant_for_uniform_measure"])
        self.assertFalse(alpha["breaks_symmetry_by_itself"])
        self.assertFalse(alpha["strict_selector_source_exported"])
        vortex = probe["vortex_torsion_audit"]
        self.assertTrue(vortex["legacy_or_numeric_vortex_terms_detected"])
        self.assertFalse(vortex["strict_twisted_vortex_source_exported"])
        self.assertFalse(vortex["strict_torsion_to_task3_provider_margin_bridge_exported"])
        self.assertFalse(vortex["legacy_beta_tors_role_transfer_allowed"])
        answers = probe["computational_answers"]
        self.assertFalse(answers["is_full_eom_exported"])
        self.assertFalse(answers["is_full_lagrangian_exported"])
        self.assertFalse(answers["does_4ln2_break_symmetry"])
        self.assertFalse(answers["does_current_strict_route_have_twisted_vortex_bridge"])
        self.assertTrue(answers["current_provider_margin_route_stopped"])
        task3 = probe["task3_g1_g3_update"]
        self.assertEqual(task3["G1_reduction_certainty"], "OPEN_UNCHANGED")
        self.assertEqual(task3["G2_nonlinear_trajectory_realism"], "CLOSED_FROM_P2301_UNCHANGED")
        self.assertEqual(task3["G3_operational_policy_rule"], "OPEN_UNCHANGED")
        theorem = probe["theorem_export"]
        self.assertEqual(sha256_json(theorem), probe["theorem_fingerprint_sha256"])
        g = data["gatekeeper_checks"]
        self.assertTrue(g["schematic_candidate_lagrangian_present"])
        self.assertTrue(g["full_tensor_eom_not_exported"])
        self.assertTrue(g["full_task3_lagrangian_not_exported"])
        self.assertTrue(g["p2030_p2033_tensor_gaps_loaded"])
        self.assertTrue(g["alpha_entropy_computation_matches_4ln2"])
        self.assertTrue(g["alpha_permutation_invariant_not_selector"])
        self.assertTrue(g["strict_twisted_vortex_bridge_not_exported"])
        self.assertTrue(g["legacy_torsion_role_not_transferred"])
        self.assertTrue(g["p2313_route_stop_loaded"])
        self.assertTrue(g["strict_g1_g3_not_updated"])
        self.assertTrue(g["no_selector_premise_added"])
        self.assertTrue(g["no_qw2191_discharge_claimed"])
        self.assertTrue(g["no_toe_closure_claimed"])


if __name__ == "__main__":
    unittest.main()
