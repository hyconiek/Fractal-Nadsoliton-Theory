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


class TestP2311S1261StrictSelfEnergyToPolicyMarginTheoremAuditProbe(unittest.TestCase):
    def test_packet(self) -> None:
        subprocess.run(
            [sys.executable, str(ROOT / "p2311_s1261_strict_self_energy_to_policy_margin_theorem_audit_probe.py")],
            check=True,
        )
        data = json.loads(
            (G / "p2311_s1261_strict_self_energy_to_policy_margin_theorem_audit_probe.json").read_text(encoding="utf-8")
        )
        self.assertEqual(data["schema_version"], "p2311_s1261_v1")
        self.assertEqual(data["packet_id"], "P2311")
        probe = data["strict_self_energy_to_policy_margin_theorem_audit_probe"]
        attempt = probe["attempted_theorem"]
        self.assertEqual(attempt["failure_mode"], "MISSING_ADM_LAPSE_SHEAR_TO_POLICY_MARGIN_RESPONSE_SOURCE")
        self.assertTrue(attempt["numeric_condition_passes"])
        self.assertTrue(attempt["p2310_replay_passes"])
        self.assertFalse(attempt["strict_theorem_proven"])
        obligations = {row["obligation_id"]: row for row in probe["theorem_obligations"]}
        self.assertTrue(obligations["O1_TYPED_SELF_ENERGY_DOMAIN"]["satisfied"])
        self.assertTrue(obligations["O2_SELF_ENERGY_NUMERIC_LIFT_AND_REPLAY"]["satisfied"])
        self.assertFalse(obligations["O3_ADM_LAPSE_SHEAR_ENERGY_SOURCE_THEOREM"]["satisfied"])
        self.assertFalse(obligations["O4_SIGNED_MONOTONE_POLICY_MARGIN_RESPONSE"]["satisfied"])
        self.assertTrue(obligations["O5_NO_SELECTOR_OR_CONVENTION_LAYER"]["satisfied"])
        audit = probe["adm_lapse_shear_evidence_audit"]
        self.assertTrue(audit["eh_lapse_shear_sign_present"])
        self.assertTrue(audit["gb_lapse_cancellation_present"])
        self.assertTrue(audit["non_gb_lapse_obstruction_present"])
        self.assertTrue(audit["non_gb_decomposition_obstruction_present"])
        verdict = probe["admissibility_verdict"]
        self.assertFalse(verdict["self_energy_to_policy_margin_theorem_proven"])
        self.assertFalse(verdict["current_export_set_can_supply_bridge_without_new_lapse_shear_source"])
        self.assertFalse(verdict["admissible_for_g1_g3_update"])
        task3 = probe["task3_g1_g3_update"]
        self.assertEqual(task3["G1_reduction_certainty"], "OPEN_UNCHANGED")
        self.assertEqual(task3["G2_nonlinear_trajectory_realism"], "CLOSED_FROM_P2301_UNCHANGED")
        self.assertEqual(task3["G3_operational_policy_rule"], "OPEN_UNCHANGED")
        self.assertEqual(sha256_json(probe["theorem_export"]), probe["theorem_fingerprint_sha256"])
        self.assertIn("O3_ADM_LAPSE_SHEAR_ENERGY_SOURCE_THEOREM", probe["theorem_export"]["proof_bits"]["failed_obligations"])
        self.assertIn("O4_SIGNED_MONOTONE_POLICY_MARGIN_RESPONSE", probe["theorem_export"]["proof_bits"]["failed_obligations"])
        g = data["gatekeeper_checks"]
        self.assertTrue(g["p1980_eh_lapse_witness_loaded"])
        self.assertTrue(g["p1981_p1983_curvature_squared_obligations_loaded"])
        self.assertTrue(g["p1984_gb_lapse_cancellation_loaded"])
        self.assertTrue(g["p1985_p1986_non_gb_obstruction_loaded"])
        self.assertTrue(g["p2300_coefficients_loaded"])
        self.assertTrue(g["p2310_self_energy_replay_loaded"])
        self.assertTrue(g["self_energy_numeric_condition_passes"])
        self.assertTrue(g["adm_lapse_shear_energy_to_policy_margin_theorem_not_exported"])
        self.assertTrue(g["signed_monotone_policy_response_not_exported"])
        self.assertTrue(g["current_export_set_blocks_self_energy_bridge_use"])
        self.assertTrue(g["strict_g1_g3_not_updated"])
        self.assertTrue(g["no_selector_premise_added"])
        self.assertTrue(g["no_qw2191_discharge_claimed"])
        self.assertTrue(g["no_legacy_kernel_role_transfer"])
        self.assertTrue(g["no_toe_closure_claimed"])


if __name__ == "__main__":
    unittest.main()
