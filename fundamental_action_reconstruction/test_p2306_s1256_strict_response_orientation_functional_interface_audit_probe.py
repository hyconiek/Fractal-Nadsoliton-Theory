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


class TestP2306S1256StrictResponseOrientationFunctionalInterfaceAuditProbe(unittest.TestCase):
    def test_packet(self) -> None:
        subprocess.run(
            [sys.executable, str(ROOT / "p2306_s1256_strict_response_orientation_functional_interface_audit_probe.py")],
            check=True,
        )
        data = json.loads(
            (G / "p2306_s1256_strict_response_orientation_functional_interface_audit_probe.json").read_text(encoding="utf-8")
        )
        self.assertEqual(data["schema_version"], "p2306_s1256_v1")
        probe = data["strict_response_orientation_functional_interface_audit_probe"]
        self.assertEqual(len(probe["required_interface_conditions"]), 4)
        self.assertTrue(all(row["current_status"].startswith(("MISSING", "NOT_SATISFIED")) for row in probe["required_interface_conditions"]))
        self.assertTrue(all(row["loaded"] for row in probe["orientation_source_audit"]))
        self.assertFalse(any(row["exports_provider_margin_response"] for row in probe["orientation_source_audit"]))
        candidates = {row["candidate_id"]: row for row in probe["candidate_response_functionals"]}
        self.assertFalse(any(row["passes_required_lift"] for row in candidates.values()))
        self.assertTrue(candidates["SIGNED_TOTAL_CANONICAL_RESPONSE"]["well_typed_on_p2300_basis"])
        self.assertLess(candidates["SIGNED_TOTAL_CANONICAL_RESPONSE"]["signed_margin_value"], 0.0)
        verdict = probe["strict_interface_verdict"]
        self.assertFalse(verdict["response_orientation_functional_exported"])
        self.assertTrue(verdict["current_export_set_refutes_closure_use"])
        self.assertEqual(probe["task3_g1_g3_update"]["G1_reduction_certainty"], "OPEN_UNCHANGED")
        self.assertEqual(probe["task3_g1_g3_update"]["G3_operational_policy_rule"], "OPEN_UNCHANGED")
        self.assertEqual(sha256_json(probe["theorem_export"]), probe["theorem_fingerprint_sha256"])
        g = data["gatekeeper_checks"]
        self.assertTrue(g["alpha_geo_strict_source_loaded"])
        self.assertTrue(g["alpha_geo_is_four_ln2_not_legacy_import"])
        self.assertTrue(g["p2300_coefficients_loaded"])
        self.assertTrue(g["p2305_underdetermination_loaded"])
        self.assertTrue(g["orientation_artifacts_audited"])
        self.assertTrue(g["no_orientation_source_exports_provider_margin_response"])
        self.assertTrue(g["candidate_response_functionals_do_not_close_required_lift"])
        self.assertTrue(g["strict_response_orientation_bridge_not_proven"])
        self.assertTrue(g["strict_g1_g3_not_updated"])
        self.assertTrue(g["no_qw2191_discharge_claimed"])
        self.assertTrue(g["no_selector_closure_claimed"])
        self.assertTrue(g["no_legacy_kernel_role_transfer"])
        self.assertTrue(g["no_toe_closure_claimed"])


if __name__ == "__main__":
    unittest.main()
