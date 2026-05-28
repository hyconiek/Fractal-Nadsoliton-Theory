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


class TestP2317S1267StrictFourierDegeneracyExistingLaneAuditProbe(unittest.TestCase):
    def test_packet(self) -> None:
        subprocess.run(
            [sys.executable, str(ROOT / "p2317_s1267_strict_fourier_degeneracy_existing_lane_audit_probe.py")],
            check=True,
        )
        data = json.loads(
            (G / "p2317_s1267_strict_fourier_degeneracy_existing_lane_audit_probe.json").read_text(
                encoding="utf-8"
            )
        )
        self.assertEqual(data["schema_version"], "p2317_s1267_v1")
        self.assertEqual(data["packet_id"], "P2317")
        self.assertEqual(
            data["result_kind"],
            "STRICT_FOURIER_PAIR_DEGENERACY_RECOMPUTATION_AND_EXISTING_SELECTOR_LANE_AUDIT_NO_G1_G3_UPDATE",
        )
        probe = data["strict_fourier_degeneracy_existing_lane_audit_probe"]
        grep_audit = probe["existing_evidence_grep_audit"]
        self.assertGreaterEqual(grep_audit["hit_count"], 10)
        self.assertGreaterEqual(len(grep_audit["top_hits"]), 5)
        kernel = probe["kernel_alone_pair_block_computation"]
        self.assertTrue(kernel["all_pair_blocks_scalar_identity"])
        self.assertTrue(kernel["p2315_pair_degeneracy_verified"])
        self.assertLess(kernel["max_scalar_identity_residual"], 1e-10)
        self.assertEqual(len(kernel["real_fourier_pair_blocks"]), 5)
        self.assertTrue(all(row["kernel_alone_degenerate"] for row in kernel["real_fourier_pair_blocks"]))
        diagonal = probe["diagonal_local_lane_o2_to_z2_computation"]
        self.assertTrue(diagonal["all_pairs_cut"])
        self.assertGreater(diagonal["min_defect_abs"], 1e-10)
        self.assertGreater(diagonal["min_eigenvalue_split"], 0.0)
        self.assertIn("not a Task-3", diagonal["lane_limitation"])
        shannon = probe["shannon_element_order_lane_o2_to_z2_computation"]
        self.assertTrue(shannon["all_defects_nonzero"])
        self.assertGreater(shannon["min_defect_abs"], 1e-10)
        self.assertGreater(shannon["min_eigenvalue_split"], 0.0)
        self.assertIn("not translation-invariant", " ".join(shannon["r_ord_translation_invariance_note"]).lower())
        comparison = probe["lane_comparison_and_route_decision"]
        self.assertFalse(comparison["kernel_alone"]["orientation_selected"])
        self.assertTrue(comparison["diagonal_local_lane"]["orientation_selected_up_to_Z2_in_lane"])
        self.assertTrue(comparison["shannon_element_order_lane"]["orientation_selected_up_to_Z2_in_lane"])
        self.assertTrue(comparison["diagonal_local_lane"]["not_a_task3_bridge"])
        self.assertTrue(comparison["shannon_element_order_lane"]["not_a_task3_bridge"])
        theorem = probe["theorem_export"]
        self.assertEqual(sha256_json(theorem), probe["theorem_fingerprint_sha256"])
        self.assertIn("already computed", theorem["claim"])
        g = data["gatekeeper_checks"]
        self.assertTrue(g["far_grep_found_existing_fourier_selector_evidence"])
        self.assertTrue(g["p2315_loaded"])
        self.assertTrue(g["p2315_pair_degeneracy_verified"])
        self.assertTrue(g["direct_kernel_pair_blocks_scalar_identity"])
        self.assertTrue(g["direct_kernel_pair_block_residual_small"])
        self.assertTrue(g["diagonal_local_lane_all_pairs_cut"])
        self.assertTrue(g["diagonal_local_lane_nonzero_min_defect"])
        self.assertTrue(g["shannon_lane_all_defects_nonzero"])
        self.assertTrue(g["shannon_lane_nonzero_min_defect"])
        self.assertTrue(g["existing_lane_exports_not_promoted_to_task3_bridge"])
        self.assertTrue(g["strict_g1_g3_not_updated"])
        self.assertTrue(g["no_selector_closure_claimed"])
        self.assertTrue(g["no_qw2191_global_discharge_claimed"])
        self.assertTrue(g["no_toe_closure_claimed"])


if __name__ == "__main__":
    unittest.main()
