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


class TestP2307S1257StrictDynamicalMarginResponseFunctionalNonderivationProbe(unittest.TestCase):
    def test_packet(self) -> None:
        subprocess.run(
            [sys.executable, str(ROOT / "p2307_s1257_strict_dynamical_margin_response_functional_nonderivation_probe.py")],
            check=True,
        )
        data = json.loads(
            (G / "p2307_s1257_strict_dynamical_margin_response_functional_nonderivation_probe.json").read_text(encoding="utf-8")
        )
        self.assertEqual(data["schema_version"], "p2307_s1257_v1")
        probe = data["strict_dynamical_margin_response_functional_nonderivation_probe"]
        audit = probe["dynamical_interface_audit"]
        self.assertEqual(audit["d_update_d_provider_lift"], "1")
        self.assertEqual(audit["d_update_d_p2300_coefficients_without_bridge"], ["0"] * 10)
        self.assertEqual(audit["chain_rule_if_lambda_equals_w_dot_c"]["d_update_d_coefficients"], [f"w{i}" for i in range(10)])
        candidates = {row["candidate_id"]: row for row in probe["candidate_functionals"]}
        self.assertTrue(candidates["DYNAMICS_NATIVE_PROVIDER_LIFT_PARAMETER"]["passes_required_lift"])
        self.assertFalse(candidates["DYNAMICS_NATIVE_PROVIDER_LIFT_PARAMETER"]["typed_on_p2300_coefficients"])
        self.assertLess(candidates["CANONICAL_SIGNED_TOTAL_FROM_P2300"]["value"], 0.0)
        self.assertFalse(any(row["admissible_as_p2300_response_functional"] for row in candidates.values()))
        verdict = probe["strict_nonderivation_verdict"]
        self.assertTrue(verdict["scalar_lift_sensitivity_derived"])
        self.assertFalse(verdict["p2300_to_margin_response_functional_derived"])
        self.assertTrue(verdict["current_dynamics_refutes_closure_use"])
        self.assertEqual(probe["task3_g1_g3_update"]["G1_reduction_certainty"], "OPEN_UNCHANGED")
        self.assertEqual(probe["task3_g1_g3_update"]["G3_operational_policy_rule"], "OPEN_UNCHANGED")
        self.assertEqual(sha256_json(probe["theorem_export"]), probe["theorem_fingerprint_sha256"])
        g = data["gatekeeper_checks"]
        self.assertTrue(g["alpha_geo_strict_source_loaded"])
        self.assertTrue(g["alpha_geo_is_four_ln2_not_legacy_import"])
        self.assertTrue(g["p2300_coefficients_loaded"])
        self.assertTrue(g["p2302_required_lift_loaded"])
        self.assertTrue(g["p2306_interface_audit_loaded"])
        self.assertTrue(g["dynamics_derivative_to_scalar_lift_is_one"])
        self.assertTrue(g["dynamics_derivative_to_coefficients_without_bridge_is_zero"])
        self.assertTrue(g["chain_rule_requires_unexported_weights"])
        self.assertTrue(g["no_candidate_derives_admissible_response_functional"])
        self.assertTrue(g["strict_response_functional_not_derived"])
        self.assertTrue(g["strict_g1_g3_not_updated"])
        self.assertTrue(g["no_qw2191_discharge_claimed"])
        self.assertTrue(g["no_selector_closure_claimed"])
        self.assertTrue(g["no_legacy_kernel_role_transfer"])
        self.assertTrue(g["no_toe_closure_claimed"])


if __name__ == "__main__":
    unittest.main()
