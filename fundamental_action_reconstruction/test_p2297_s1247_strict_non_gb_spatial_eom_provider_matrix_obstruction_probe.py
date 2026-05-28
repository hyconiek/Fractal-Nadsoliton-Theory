from __future__ import annotations
import json
import subprocess
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent
G = ROOT / "generated"


class TestP2297S1247StrictNonGbSpatialEomProviderMatrixObstructionProbe(unittest.TestCase):
    def test_packet(self) -> None:
        subprocess.run(
            [sys.executable, str(ROOT / "p2297_s1247_strict_non_gb_spatial_eom_provider_matrix_obstruction_probe.py")],
            check=True,
        )
        data = json.loads(
            (G / "p2297_s1247_strict_non_gb_spatial_eom_provider_matrix_obstruction_probe.json").read_text(encoding="utf-8")
        )
        self.assertEqual(data["schema_version"], "p2297_s1247_v1")
        probe = data["strict_non_gb_spatial_eom_provider_matrix_obstruction_probe"]
        results = probe["provider_matrix_results"]
        self.assertFalse(results["strict_core_minimal_provider"]["matrix_report"]["consistent"])
        self.assertFalse(results["p1990_augmented_provider_non_strict"]["matrix_report"]["consistent"])
        self.assertTrue(results["formal_full_residual_basis_provider"]["exact_reconstruction_zero"])
        self.assertIn("NOT_A_LEGAL_STRICT_PROVIDER", results["formal_full_residual_basis_provider"]["admissibility_verdict"])
        g = data["gatekeeper_checks"]
        self.assertTrue(g["tracefree_spatial_component_sum_zero"])
        self.assertTrue(g["strict_core_matrix_inconsistent"])
        self.assertTrue(g["augmented_non_strict_matrix_inconsistent"])
        self.assertTrue(g["formal_full_basis_reconstructs_zero"])
        self.assertTrue(g["formal_full_basis_marked_non_admissible"])
        self.assertTrue(g["p1985_obstruction_preserved"])
        self.assertTrue(g["no_legacy_transfer_used"])
        self.assertTrue(g["qw2191_not_discharged"])
        self.assertTrue(g["no_selector_closure_claimed"])
        self.assertTrue(g["no_toe_closure_claimed"])


if __name__ == "__main__":
    unittest.main()
