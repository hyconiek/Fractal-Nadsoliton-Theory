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


class TestP2310S1260StrictSelfEnergyResponseSourceAuditAndReplayProbe(unittest.TestCase):
    def test_packet(self) -> None:
        subprocess.run(
            [sys.executable, str(ROOT / "p2310_s1260_strict_self_energy_response_source_audit_and_replay_probe.py")],
            check=True,
        )
        data = json.loads(
            (G / "p2310_s1260_strict_self_energy_response_source_audit_and_replay_probe.json").read_text(encoding="utf-8")
        )
        self.assertEqual(data["schema_version"], "p2310_s1260_v1")
        probe = data["strict_self_energy_response_source_audit_and_replay_probe"]
        attempt = probe["self_energy_variational_attempt"]
        self.assertFalse(attempt["target_calibrated"])
        self.assertFalse(attempt["strict_admissible"])
        self.assertGreater(attempt["lambda_numeric"], 0.0068)
        replay = probe["self_energy_replay"]["summary"]
        self.assertTrue(replay["all_rows_meet_target"])
        self.assertGreater(replay["worst_margin_to_target"], 0.0)
        candidates = {row["candidate_id"]: row for row in probe["source_candidates"]}
        self.assertFalse(candidates["ADM_COEFFICIENT_SELF_ENERGY_NORM_SQUARED"]["target_calibrated"])
        self.assertTrue(candidates["ADM_COEFFICIENT_SELF_ENERGY_NORM_SQUARED"]["replay_all_rows_meet_target"])
        self.assertFalse(any(row["strict_admissible"] for row in candidates.values()))
        audit = probe["adm_source_audit"]
        self.assertTrue(audit["p1979_no_exported_lapse_shear_provider_certificate"])
        verdict = probe["strict_admissibility_verdict"]
        self.assertTrue(verdict["non_target_lambda_candidate_found"])
        self.assertTrue(verdict["non_target_lambda_candidate_replay_passes"])
        self.assertFalse(verdict["adm_bianchi_margin_source_exported"])
        self.assertFalse(verdict["admissible_for_g1_g3_update"])
        self.assertEqual(probe["task3_g1_g3_update"]["G1_reduction_certainty"], "OPEN_UNCHANGED")
        self.assertEqual(probe["task3_g1_g3_update"]["G3_operational_policy_rule"], "OPEN_UNCHANGED")
        self.assertEqual(sha256_json(probe["theorem_export"]), probe["theorem_fingerprint_sha256"])
        g = data["gatekeeper_checks"]
        self.assertTrue(g["alpha_geo_strict_source_loaded"])
        self.assertTrue(g["alpha_geo_is_four_ln2_not_legacy_import"])
        self.assertTrue(g["p2300_coefficients_loaded"])
        self.assertTrue(g["p2302_required_lift_loaded"])
        self.assertTrue(g["p2309_quarantine_loaded"])
        self.assertTrue(g["self_energy_candidate_not_target_calibrated"])
        self.assertTrue(g["self_energy_candidate_passes_replay"])
        self.assertTrue(g["p1979_blocks_adm_lapse_shear_source_export"])
        self.assertTrue(g["no_source_candidate_strict_admissible"])
        self.assertTrue(g["strict_self_energy_bridge_not_claimed"])
        self.assertTrue(g["strict_g1_g3_not_updated"])
        self.assertTrue(g["no_qw2191_discharge_claimed"])
        self.assertTrue(g["no_selector_closure_claimed"])
        self.assertTrue(g["no_legacy_kernel_role_transfer"])
        self.assertTrue(g["no_toe_closure_claimed"])


if __name__ == "__main__":
    unittest.main()
