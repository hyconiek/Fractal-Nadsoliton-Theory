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


class TestP2324S1274StrictAxisBranchSusceptibilityNonpromotionProbe(unittest.TestCase):
    def test_packet(self) -> None:
        subprocess.run(
            [sys.executable, str(ROOT / "p2324_s1274_strict_axis_branch_susceptibility_nonpromotion_probe.py")],
            check=True,
        )
        data = json.loads(
            (G / "p2324_s1274_strict_axis_branch_susceptibility_nonpromotion_probe.json").read_text(
                encoding="utf-8"
            )
        )
        self.assertEqual(data["schema_version"], "p2324_s1274_v1")
        self.assertEqual(data["packet_id"], "P2324")
        self.assertEqual(data["result_kind"], "STRICT_AXIS_BRANCH_SUSCEPTIBILITY_NONPROMOTION_AUDIT_NO_G1_G3_UPDATE")
        probe = data["strict_axis_branch_susceptibility_nonpromotion_probe"]
        self.assertGreaterEqual(probe["far_susceptibility_grep_audit"]["hit_count"], 5)
        rows = probe["susceptibility_rows"]
        self.assertEqual(len(rows), 5)
        for row in rows:
            self.assertGreater(row["spread_per_unit_mu"], 0.0)
            self.assertTrue(row["splits_degenerate_branch_orbit_if_mu_given"])
            self.assertFalse(row["strict_internal_source_for_mu_exported"])
            self.assertAlmostEqual(row["required_lift_spread_check"], 0.0068)
            self.assertGreater(row["mu_needed_for_required_lift_spread"], 0.0)
            self.assertIn(0, row["selected_branch_indices_for_positive_mu"])
        cert = probe["susceptibility_certificate"]
        self.assertTrue(cert["all_pairs_split_by_template_if_mu_given"])
        self.assertFalse(cert["strict_internal_mu_source_exported"])
        self.assertLess(cert["max_required_lift_spread_residual"], 1e-12)
        self.assertEqual(len(cert["lay_symptom_summary"]), 3)
        update = probe["bridge_obligation_update"]
        self.assertAlmostEqual(update["required_lift_per_step"], 0.0068)
        self.assertTrue(update["signed_perturbation_susceptibility_exported"])
        self.assertFalse(update["strict_internal_signed_perturbation_exported"])
        self.assertFalse(update["susceptibility_fills_any_missing_p2318_field"])
        self.assertEqual(len(update["fields_still_unfilled_after_susceptibility_trace"]), 6)
        theorem = probe["theorem_export"]
        self.assertEqual(sha256_json(theorem), probe["theorem_fingerprint_sha256"])
        self.assertIn("high susceptibility", theorem["formal_statement"])
        g = data["gatekeeper_checks"]
        self.assertTrue(g["far_susceptibility_grep_hits_found"])
        self.assertTrue(g["p2323_loaded"])
        self.assertTrue(g["p2318_loaded"])
        self.assertTrue(g["all_pairs_split_by_template_if_mu_given"])
        self.assertTrue(g["required_lift_spread_reconstructed"])
        self.assertTrue(g["strict_internal_mu_source_not_exported"])
        self.assertTrue(g["susceptibility_not_promoted_to_bridge"])
        self.assertTrue(g["p2318_bridge_fields_still_unfilled"])
        self.assertTrue(g["strict_g1_g3_not_updated"])
        self.assertTrue(g["no_selector_closure_claimed"])
        self.assertTrue(g["no_qw2191_discharge_claimed"])
        self.assertTrue(g["no_legacy_kernel_role_transfer"])
        self.assertTrue(g["no_toe_closure_claimed"])


if __name__ == "__main__":
    unittest.main()
