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


class TestP2313S1263StrictProviderMarginRouteStopAndSourceClassContractProbe(unittest.TestCase):
    def test_packet(self) -> None:
        subprocess.run(
            [sys.executable, str(ROOT / "p2313_s1263_strict_provider_margin_route_stop_and_source_class_contract_probe.py")],
            check=True,
        )
        data = json.loads(
            (G / "p2313_s1263_strict_provider_margin_route_stop_and_source_class_contract_probe.json").read_text(encoding="utf-8")
        )
        self.assertEqual(data["schema_version"], "p2313_s1263_v1")
        self.assertEqual(data["packet_id"], "P2313")
        probe = data["strict_provider_margin_route_stop_and_source_class_contract_probe"]
        blockers = probe["blocker_stack"]
        self.assertEqual([row["packet"] for row in blockers], ["P2308", "P2309", "P2310", "P2311", "P2312"])
        self.assertTrue(all(row["active"] for row in blockers))
        audit = probe["dynamical_orientation_constraint_audit"]
        self.assertFalse(audit["constraint_exported_now"])
        self.assertTrue(audit["positive_negative_witness_available"])
        self.assertFalse(any(row["strict_admissible"] for row in audit["candidate_constraints"]))
        contract = probe["new_source_class_admission_contract"]
        self.assertEqual([row["criterion_id"] for row in contract], ["NSC1_TYPED_DOMAIN", "NSC2_INTERNAL_STRICT_ORIGIN", "NSC3_SIGNED_MONOTONE_RESPONSE", "NSC4_SELECTOR_FREE_OR_EXPLICITLY_NONSTRICT", "NSC5_REPLAY_AND_GATE_ROWS", "NSC6_GUARDRAILS"])
        self.assertTrue(all(row["mandatory"] for row in contract))
        candidates = probe["candidate_next_source_classes"]
        self.assertFalse(any(row["admitted_now"] for row in candidates))
        verdict = probe["route_stop_verdict"]
        self.assertTrue(verdict["current_provider_margin_bridge_route_stopped"])
        self.assertTrue(verdict["not_universal_no_go"])
        self.assertTrue(verdict["new_source_class_required"])
        self.assertFalse(verdict["admissible_for_g1_g3_update"])
        task3 = probe["task3_g1_g3_update"]
        self.assertEqual(task3["G1_reduction_certainty"], "OPEN_UNCHANGED")
        self.assertEqual(task3["G2_nonlinear_trajectory_realism"], "CLOSED_FROM_P2301_UNCHANGED")
        self.assertEqual(task3["G3_operational_policy_rule"], "OPEN_UNCHANGED")
        theorem = probe["theorem_export"]
        self.assertTrue(theorem["proof_bits"]["all_blockers_active"])
        self.assertTrue(theorem["proof_bits"]["new_source_class_required"])
        self.assertEqual(sha256_json(theorem), probe["theorem_fingerprint_sha256"])
        g = data["gatekeeper_checks"]
        self.assertTrue(g["p2308_current_interface_blocker_loaded"])
        self.assertTrue(g["p2309_target_calibration_quarantine_loaded"])
        self.assertTrue(g["p2310_self_energy_quarantine_loaded"])
        self.assertTrue(g["p2311_self_energy_theorem_blocker_loaded"])
        self.assertTrue(g["p2312_sign_indefinite_source_blocker_loaded"])
        self.assertTrue(g["all_current_route_blockers_active"])
        self.assertTrue(g["no_current_dynamical_orientation_constraint_exported"])
        self.assertTrue(g["new_source_class_contract_exported"])
        self.assertTrue(g["current_route_stop_scope_not_universal_no_go"])
        self.assertTrue(g["strict_g1_g3_not_updated"])
        self.assertTrue(g["no_selector_premise_added"])
        self.assertTrue(g["no_qw2191_discharge_claimed"])
        self.assertTrue(g["no_legacy_kernel_role_transfer"])
        self.assertTrue(g["no_toe_closure_claimed"])


if __name__ == "__main__":
    unittest.main()
