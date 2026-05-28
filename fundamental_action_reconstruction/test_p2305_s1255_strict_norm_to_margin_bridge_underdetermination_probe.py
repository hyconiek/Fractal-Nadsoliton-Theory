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


class TestP2305S1255StrictNormToMarginBridgeUnderdeterminationProbe(unittest.TestCase):
    def test_packet(self) -> None:
        subprocess.run(
            [sys.executable, str(ROOT / "p2305_s1255_strict_norm_to_margin_bridge_underdetermination_probe.py")],
            check=True,
        )
        data = json.loads(
            (G / "p2305_s1255_strict_norm_to_margin_bridge_underdetermination_probe.json").read_text(encoding="utf-8")
        )
        self.assertEqual(data["schema_version"], "p2305_s1255_v1")
        probe = data["strict_norm_to_margin_bridge_underdetermination_probe"]
        norms = probe["norm_summary"]
        self.assertGreater(norms["l1_norm"], norms["required_lift"])
        self.assertGreater(norms["l2_norm"], norms["required_lift"])
        self.assertLess(norms["signed_total"], 0.0)
        witnesses = {row["witness_id"]: row for row in probe["response_functional_underdetermination_witnesses"]}
        self.assertGreater(witnesses["UNIT_ALIGNED_POSITIVE_RESPONSE"]["lift_value"], norms["required_lift"])
        self.assertLess(witnesses["UNIT_ANTI_ALIGNED_NEGATIVE_RESPONSE"]["lift_value"], 0.0)
        self.assertAlmostEqual(witnesses["UNIT_ORTHOGONAL_ZERO_RESPONSE"]["lift_value"], 0.0, places=12)
        verdict = probe["strict_bridge_verdict"]
        self.assertFalse(verdict["norm_to_margin_bridge_proven"])
        self.assertTrue(verdict["norm_to_margin_bridge_refuted_as_norm_only_theorem"])
        self.assertEqual(probe["task3_g1_g3_update"]["G1_reduction_certainty"], "OPEN_UNCHANGED")
        self.assertEqual(probe["task3_g1_g3_update"]["G3_operational_policy_rule"], "OPEN_UNCHANGED")
        self.assertEqual(sha256_json(probe["theorem_export"]), probe["theorem_fingerprint_sha256"])
        g = data["gatekeeper_checks"]
        self.assertTrue(g["alpha_geo_strict_source_loaded"])
        self.assertTrue(g["alpha_geo_is_four_ln2_not_legacy_import"])
        self.assertTrue(g["p2300_coefficients_loaded"])
        self.assertTrue(g["p2302_required_lift_loaded"])
        self.assertTrue(g["p2304_direct_floor_refutation_loaded"])
        self.assertTrue(g["norms_numerically_exceed_required_lift"])
        self.assertTrue(g["anti_aligned_witness_refutes_norm_only_implication"])
        self.assertTrue(g["orthogonal_witness_refutes_norm_only_implication"])
        self.assertTrue(g["strict_norm_to_margin_bridge_not_proven"])
        self.assertTrue(g["strict_g1_g3_not_updated"])
        self.assertTrue(g["no_qw2191_discharge_claimed"])
        self.assertTrue(g["no_selector_closure_claimed"])
        self.assertTrue(g["no_legacy_kernel_role_transfer"])
        self.assertTrue(g["no_toe_closure_claimed"])


if __name__ == "__main__":
    unittest.main()
