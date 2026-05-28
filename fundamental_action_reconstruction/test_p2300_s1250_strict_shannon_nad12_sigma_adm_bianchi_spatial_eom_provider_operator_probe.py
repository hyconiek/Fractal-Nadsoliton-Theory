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


class TestP2300S1250StrictShannonNad12SigmaAdmBianchiSpatialEomProviderOperatorProbe(unittest.TestCase):
    def test_packet(self) -> None:
        subprocess.run(
            [sys.executable, str(ROOT / "p2300_s1250_strict_shannon_nad12_sigma_adm_bianchi_spatial_eom_provider_operator_probe.py")],
            check=True,
        )
        data = json.loads(
            (G / "p2300_s1250_strict_shannon_nad12_sigma_adm_bianchi_spatial_eom_provider_operator_probe.json").read_text(encoding="utf-8")
        )
        self.assertEqual(data["schema_version"], "p2300_s1250_v1")
        probe = data["strict_shannon_nad12_sigma_adm_bianchi_spatial_eom_provider_operator_probe"]
        self.assertEqual(probe["operator_basis_policy"]["column_count"], 10)
        self.assertFalse(probe["operator_basis_policy"]["uses_formal_full_residual_basis"])
        self.assertFalse(probe["operator_basis_policy"]["uses_explicit_selector_premise"])
        report = probe["provider_matrix_report"]
        self.assertTrue(report["consistent"])
        self.assertEqual(report["rank_A"], report["rank_augmented"])
        self.assertTrue(probe["solution_space"]["exact_reconstruction_zero"])
        task3 = probe["task3_g1_g2_g3_impact"]
        self.assertEqual(task3["gap_status_before_p2300"], task3["gap_status_after_p2300"])
        self.assertTrue(all(status == "OPEN" for status in task3["gap_status_after_p2300"].values()))
        self.assertEqual(sha256_json(probe["theorem_export"]), probe["theorem_fingerprint_sha256"])
        g = data["gatekeeper_checks"]
        self.assertTrue(g["alpha_geo_strict_source_loaded"])
        self.assertTrue(g["alpha_geo_is_four_ln2_not_legacy_import"])
        self.assertTrue(g["p1988_family_scope_loaded"])
        self.assertTrue(g["p2299_no_column_blocker_was_real"])
        self.assertTrue(g["new_spatial_eom_operator_columns_exported"])
        self.assertTrue(g["operator_basis_not_formal_residual_copy"])
        self.assertTrue(g["provider_matrix_consistent"])
        self.assertTrue(g["exact_reconstruction_zero"])
        self.assertTrue(g["canonical_solution_exported"])
        self.assertTrue(g["strict_selector_premise_not_used"])
        self.assertTrue(g["task3_g1_g2_g3_not_closed_by_p2300"])
        self.assertTrue(g["no_qw2191_discharge_claimed"])
        self.assertTrue(g["no_legacy_kernel_role_transfer"])
        self.assertTrue(g["no_toe_closure_claimed"])


if __name__ == "__main__":
    unittest.main()
