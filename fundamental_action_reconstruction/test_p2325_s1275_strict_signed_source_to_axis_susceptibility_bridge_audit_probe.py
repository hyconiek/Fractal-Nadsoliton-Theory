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


class TestP2325S1275StrictSignedSourceToAxisSusceptibilityBridgeAuditProbe(unittest.TestCase):
    def test_packet(self) -> None:
        subprocess.run(
            [sys.executable, str(ROOT / "p2325_s1275_strict_signed_source_to_axis_susceptibility_bridge_audit_probe.py")],
            check=True,
        )
        data = json.loads(
            (G / "p2325_s1275_strict_signed_source_to_axis_susceptibility_bridge_audit_probe.json").read_text(
                encoding="utf-8"
            )
        )
        self.assertEqual(data["schema_version"], "p2325_s1275_v1")
        self.assertEqual(data["packet_id"], "P2325")
        self.assertEqual(data["result_kind"], "STRICT_SIGNED_SOURCE_TO_AXIS_SUSCEPTIBILITY_BRIDGE_AUDIT_NO_G1_G3_UPDATE")
        probe = data["strict_signed_source_to_axis_susceptibility_bridge_audit_probe"]
        self.assertGreaterEqual(probe["far_signed_source_grep_audit"]["hit_count"], 5)
        self.assertEqual(len(probe["required_signed_source_fields"]), 6)
        self.assertEqual(len(probe["candidate_source_rows"]), 5)
        matrix = probe["predicate_matrix"]
        self.assertEqual(len(matrix), 5)
        self.assertTrue(all(not row["admissible_for_p2324_signed_mu_bridge"] for row in matrix))
        self.assertTrue(all(row["failed_predicates"] for row in matrix))
        p731 = next(row for row in probe["candidate_source_rows"] if row["candidate_id"] == "P731_W_BREAK_PAIR12_ANTISYMMETRIC_SCORES")
        self.assertTrue(p731["theta_ref_or_branch_reference"])
        self.assertFalse(p731["strict_internal_signed_mu_source"])
        cert = probe["bridge_audit_certificate"]
        self.assertEqual(cert["candidate_count"], 5)
        self.assertEqual(cert["admissible_candidate_count"], 0)
        self.assertEqual(cert["closest_candidate"], "P731_W_BREAK_PAIR12_ANTISYMMETRIC_SCORES")
        self.assertTrue(cert["all_candidates_fail_at_least_one_required_predicate"])
        self.assertEqual(len(cert["p2324_required_mu_by_pair_loaded"]), 5)
        update = probe["bridge_obligation_update"]
        self.assertFalse(update["signed_source_to_mu_bridge_exported"])
        self.assertEqual(len(update["fields_still_unfilled_after_signed_source_audit"]), 6)
        self.assertFalse(update["g1_g3_update_allowed"])
        theorem = probe["theorem_export"]
        self.assertEqual(sha256_json(theorem), probe["theorem_fingerprint_sha256"])
        self.assertIn("no admissible signed-source bridge", theorem["formal_statement"])
        g = data["gatekeeper_checks"]
        self.assertTrue(g["far_signed_source_grep_hits_found"])
        self.assertTrue(g["p2324_loaded"])
        self.assertTrue(g["p2318_loaded"])
        self.assertTrue(g["p1967_loaded"])
        self.assertTrue(g["candidate_class_nonempty"])
        self.assertTrue(g["p2324_mu_values_loaded"])
        self.assertTrue(g["no_admissible_signed_source_bridge"])
        self.assertTrue(g["all_candidates_fail_at_least_one_required_predicate"])
        self.assertTrue(g["p731_closest_candidate_not_promoted"])
        self.assertTrue(g["p2318_bridge_fields_still_unfilled"])
        self.assertTrue(g["strict_g1_g3_not_updated"])
        self.assertTrue(g["no_selector_closure_claimed"])
        self.assertTrue(g["no_qw2191_discharge_claimed"])
        self.assertTrue(g["no_legacy_kernel_role_transfer"])
        self.assertTrue(g["no_toe_closure_claimed"])


if __name__ == "__main__":
    unittest.main()
