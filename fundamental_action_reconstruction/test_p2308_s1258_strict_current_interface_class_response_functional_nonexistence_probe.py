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


class TestP2308S1258StrictCurrentInterfaceClassResponseFunctionalNonexistenceProbe(unittest.TestCase):
    def test_packet(self) -> None:
        subprocess.run(
            [sys.executable, str(ROOT / "p2308_s1258_strict_current_interface_class_response_functional_nonexistence_probe.py")],
            check=True,
        )
        data = json.loads(
            (G / "p2308_s1258_strict_current_interface_class_response_functional_nonexistence_probe.json").read_text(encoding="utf-8")
        )
        self.assertEqual(data["schema_version"], "p2308_s1258_v1")
        probe = data["strict_current_interface_class_response_functional_nonexistence_probe"]
        cls = probe["current_interface_class"]
        self.assertGreaterEqual(len(cls), 8)
        self.assertTrue(all(row["failed_predicates"] for row in cls))
        self.assertFalse(any(row["passes_all_required_predicates"] for row in cls))
        cert = probe["class_nonexistence_certificate"]
        self.assertEqual(cert["candidate_count"], len(cls))
        self.assertEqual(cert["admissible_candidate_count"], 0)
        self.assertTrue(cert["nonexistence_over_current_class"])
        verdict = probe["strict_nonexistence_verdict"]
        self.assertFalse(verdict["lambda_equals_R_of_c_exported_in_current_class"])
        self.assertTrue(verdict["nonexistence_proven_for_current_interface_class"])
        self.assertFalse(verdict["universal_future_nonexistence_claimed"])
        self.assertEqual(probe["task3_g1_g3_update"]["G1_reduction_certainty"], "OPEN_UNCHANGED")
        self.assertEqual(probe["task3_g1_g3_update"]["G3_operational_policy_rule"], "OPEN_UNCHANGED")
        self.assertEqual(sha256_json(probe["theorem_export"]), probe["theorem_fingerprint_sha256"])
        g = data["gatekeeper_checks"]
        self.assertTrue(g["alpha_geo_strict_source_loaded"])
        self.assertTrue(g["alpha_geo_is_four_ln2_not_legacy_import"])
        self.assertTrue(g["p2300_coefficients_loaded"])
        self.assertTrue(g["p2302_required_lift_loaded"])
        self.assertTrue(g["p2306_candidates_loaded"])
        self.assertTrue(g["p2307_candidates_loaded"])
        self.assertTrue(g["all_candidates_fail_at_least_one_required_predicate"])
        self.assertTrue(g["no_admissible_current_class_response_functional"])
        self.assertTrue(g["theorem_scoped_to_current_export_class"])
        self.assertTrue(g["strict_response_functional_not_exported"])
        self.assertTrue(g["strict_g1_g3_not_updated"])
        self.assertTrue(g["no_qw2191_discharge_claimed"])
        self.assertTrue(g["no_selector_closure_claimed"])
        self.assertTrue(g["no_legacy_kernel_role_transfer"])
        self.assertTrue(g["no_toe_closure_claimed"])


if __name__ == "__main__":
    unittest.main()
