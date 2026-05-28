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


class TestP2303S1253StrictProviderToMarginBridgeBoundAuditProbe(unittest.TestCase):
    def test_packet(self) -> None:
        subprocess.run(
            [sys.executable, str(ROOT / "p2303_s1253_strict_provider_to_margin_bridge_bound_audit_probe.py")],
            check=True,
        )
        data = json.loads(
            (G / "p2303_s1253_strict_provider_to_margin_bridge_bound_audit_probe.json").read_text(encoding="utf-8")
        )
        self.assertEqual(data["schema_version"], "p2303_s1253_v1")
        probe = data["strict_provider_to_margin_bridge_bound_audit_probe"]
        self.assertFalse(probe["strict_bridge_verdict"]["provider_to_margin_bridge_proven"])
        self.assertTrue(probe["strict_bridge_verdict"]["provider_to_margin_bridge_refuted_for_current_direct_contracts"])
        required = probe["required_lift"]["value"]
        direct = next(row for row in probe["bridge_contract_audit"] if row["contract_id"] == "STRICT_SINGLE_CHANNEL_FLOOR")
        self.assertLess(direct["bound_numeric"], required)
        self.assertFalse(direct["passes_required_lift"])
        positive_sum = next(row for row in probe["bridge_contract_audit"] if row["contract_id"] == "POSITIVE_CHANNEL_SUM")
        self.assertTrue(positive_sum["passes_required_lift"])
        self.assertFalse(positive_sum["strict_admissible_without_extra_premise"])
        self.assertEqual(sha256_json(probe["theorem_export"]), probe["theorem_fingerprint_sha256"])
        g = data["gatekeeper_checks"]
        self.assertTrue(g["alpha_geo_strict_source_loaded"])
        self.assertTrue(g["alpha_geo_is_four_ln2_not_legacy_import"])
        self.assertTrue(g["p2300_coefficients_loaded"])
        self.assertTrue(g["p2302_required_lift_loaded"])
        self.assertTrue(g["direct_floor_below_required_lift"])
        self.assertTrue(g["signed_total_not_positive_lift"])
        self.assertTrue(g["strict_direct_bridge_not_proven"])
        self.assertTrue(g["strict_direct_bridge_refuted_for_current_contracts"])
        self.assertTrue(g["sufficient_norms_marked_non_admissible_without_new_bridge"])
        self.assertTrue(g["strict_g1_g3_not_updated"])
        self.assertTrue(g["no_qw2191_discharge_claimed"])
        self.assertTrue(g["no_selector_closure_claimed"])
        self.assertTrue(g["no_legacy_kernel_role_transfer"])
        self.assertTrue(g["no_toe_closure_claimed"])


if __name__ == "__main__":
    unittest.main()
