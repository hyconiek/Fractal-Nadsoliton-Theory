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


class TestP2299S1249StrictShannonProviderSourceAttemptAndNonStrictSelectorBranchProbe(unittest.TestCase):
    def test_packet(self) -> None:
        subprocess.run(
            [sys.executable, str(ROOT / "p2299_s1249_strict_shannon_provider_source_attempt_and_non_strict_selector_branch_probe.py")],
            check=True,
        )
        data = json.loads(
            (G / "p2299_s1249_strict_shannon_provider_source_attempt_and_non_strict_selector_branch_probe.json").read_text(encoding="utf-8")
        )
        self.assertEqual(data["schema_version"], "p2299_s1249_v1")
        probe = data["strict_shannon_provider_source_attempt_and_non_strict_selector_branch_probe"]
        attempt = probe["strict_shannon_provider_attempt"]
        self.assertTrue(attempt["uses_alpha_geo_strict_derived_v1"])
        self.assertEqual(attempt["full_residual_family_cancellation_status"], "FAILS_CURRENT_STRICT_PROVIDER_TEST")
        self.assertTrue(attempt["does_not_copy_residual_basis"])
        self.assertFalse(attempt["new_spatial_eom_operator_columns_exported"])
        self.assertTrue(attempt["scaled_rank_invariance_certificate"]["inconsistency_preserved"])
        branch = probe["non_strict_selector_premise_branch"]
        self.assertEqual(branch["status"], "FORMALIZED_AS_NON_STRICT_BRANCH_ONLY")
        self.assertFalse(branch["may_close_strict_G1_G2_G3"])
        task3 = probe["task3_g1_g2_g3_impact"]
        self.assertEqual(task3["gap_status_after_p2299"], task3["gap_status_before_p2299"])
        self.assertTrue(all(status == "OPEN" for status in task3["gap_status_after_p2299"].values()))
        self.assertEqual(sha256_json(probe["theorem_export"]), probe["theorem_fingerprint_sha256"])
        g = data["gatekeeper_checks"]
        self.assertTrue(g["alpha_geo_strict_source_loaded"])
        self.assertTrue(g["alpha_geo_is_four_ln2_not_legacy_import"])
        self.assertTrue(g["shannon_weights_nonuniform"])
        self.assertTrue(g["branch_B_noncyclic_candidate_seen"])
        self.assertTrue(g["candidate_does_not_copy_residual_basis"])
        self.assertTrue(g["no_new_spatial_eom_operator_columns_exported"])
        self.assertTrue(g["scaled_strict_core_inconsistency_preserved"])
        self.assertTrue(g["p2298_no_current_legal_strict_provider_confirmed"])
        self.assertTrue(g["strict_shannon_provider_attempt_fails_full_residual"])
        self.assertTrue(g["non_strict_selector_branch_formalized"])
        self.assertTrue(g["task3_g1_g2_g3_not_closed_by_p2299"])
        self.assertTrue(g["qw2191_not_discharged"])
        self.assertTrue(g["no_selector_closure_claimed"])
        self.assertTrue(g["no_legacy_kernel_role_transfer"])
        self.assertTrue(g["no_toe_closure_claimed"])


if __name__ == "__main__":
    unittest.main()
