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


class TestP2315S1265StrictSchematicLagrangianEomKernelSpectrumProbe(unittest.TestCase):
    def test_packet(self) -> None:
        subprocess.run(
            [sys.executable, str(ROOT / "p2315_s1265_strict_schematic_lagrangian_eom_kernel_spectrum_probe.py")],
            check=True,
        )
        data = json.loads(
            (G / "p2315_s1265_strict_schematic_lagrangian_eom_kernel_spectrum_probe.json").read_text(encoding="utf-8")
        )
        self.assertEqual(data["schema_version"], "p2315_s1265_v1")
        self.assertEqual(data["packet_id"], "P2315")
        probe = data["strict_schematic_lagrangian_eom_kernel_spectrum_probe"]
        eom = probe["schematic_eom_export"]
        self.assertTrue(eom["potential_term_unspecified"])
        self.assertTrue(eom["metric_tensor_variation_missing"])
        self.assertFalse(eom["full_tensor_eom_exported"])
        spectrum = probe["kernel_spectrum_audit"]
        self.assertTrue(spectrum["matrix_properties"]["real_symmetric"])
        self.assertTrue(spectrum["matrix_properties"]["circulant"])
        self.assertEqual(len(spectrum["first_row"]), 12)
        self.assertEqual(len(spectrum["fourier_eigenvalues"]), 12)
        degeneracy = spectrum["pair_degeneracy_report"]
        self.assertTrue(degeneracy["all_pair_planes_degenerate"])
        self.assertEqual(len(degeneracy["pair_planes"]), 5)
        self.assertLess(degeneracy["max_pair_difference"], 1e-12)
        alpha_degeneracy = spectrum["alpha_scaled_pair_degeneracy_report"]
        self.assertTrue(alpha_degeneracy["all_pair_planes_degenerate"])
        verdict = probe["computational_verdict"]
        self.assertTrue(verdict["schematic_eom_derived"])
        self.assertFalse(verdict["full_eom_derived"])
        self.assertFalse(verdict["full_lagrangian_derived"])
        self.assertFalse(verdict["kernel_breaks_pair_plane_degeneracy"])
        self.assertFalse(verdict["alpha_4ln2_breaks_pair_plane_degeneracy"])
        self.assertFalse(verdict["selector_source_exported"])
        self.assertFalse(verdict["admissible_for_g1_g3_update"])
        task3 = probe["task3_g1_g3_update"]
        self.assertEqual(task3["G1_reduction_certainty"], "OPEN_UNCHANGED")
        self.assertEqual(task3["G2_nonlinear_trajectory_realism"], "CLOSED_FROM_P2301_UNCHANGED")
        self.assertEqual(task3["G3_operational_policy_rule"], "OPEN_UNCHANGED")
        theorem = probe["theorem_export"]
        self.assertEqual(sha256_json(theorem), probe["theorem_fingerprint_sha256"])
        g = data["gatekeeper_checks"]
        self.assertTrue(g["release7_schematic_lagrangian_detected"])
        self.assertTrue(g["schematic_eom_derived"])
        self.assertTrue(g["kernel_matrix_real_symmetric"])
        self.assertTrue(g["kernel_matrix_circulant"])
        self.assertTrue(g["fourier_pair_degeneracy_verified"])
        self.assertTrue(g["alpha_scaling_preserves_pair_degeneracy"])
        self.assertTrue(g["full_eom_not_derived"])
        self.assertTrue(g["full_lagrangian_not_derived"])
        self.assertTrue(g["selector_source_not_exported"])
        self.assertTrue(g["p2314_inventory_loaded"])
        self.assertTrue(g["strict_g1_g3_not_updated"])
        self.assertTrue(g["no_selector_premise_added"])
        self.assertTrue(g["no_qw2191_discharge_claimed"])
        self.assertTrue(g["no_toe_closure_claimed"])


if __name__ == "__main__":
    unittest.main()
