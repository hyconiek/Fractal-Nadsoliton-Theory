from __future__ import annotations

import hashlib
import json
import subprocess
import sys
import unittest
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parent
SCRIPT = ROOT / "p2332_s1282_strict_bianchi_spatial_tensor_component_table_gb_lift_probe.py"
OUT = ROOT / "generated" / "p2332_s1282_strict_bianchi_spatial_tensor_component_table_gb_lift_probe.json"
MD = ROOT / "generated" / "p2332_s1282_strict_bianchi_spatial_tensor_component_table_gb_lift_probe.md"


def sha256_json(obj: Any) -> str:
    return hashlib.sha256(
        json.dumps(obj, sort_keys=True, ensure_ascii=False, separators=(",", ":")).encode("utf-8")
    ).hexdigest()


class P2332BianchiSpatialTensorTableGBLiftTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        subprocess.run([sys.executable, str(SCRIPT)], check=True, cwd=ROOT.parent)
        cls.payload = json.loads(OUT.read_text(encoding="utf-8"))
        cls.probe = cls.payload["strict_bianchi_spatial_tensor_component_table_gb_lift_probe"]

    def test_packet_identity_and_outputs(self) -> None:
        self.assertEqual(self.payload["schema_version"], "p2332_s1282_v1")
        self.assertEqual(self.payload["packet_id"], "P2332")
        self.assertEqual(self.payload["stage_id"], "S1282")
        self.assertTrue(OUT.exists())
        self.assertTrue(MD.exists())

    def test_component_table_and_same_basis_target(self) -> None:
        summary = self.probe["component_table_summary"]
        self.assertEqual(summary["component_entry_count"], 12)
        self.assertTrue(summary["all_tracefree_sums_zero"])
        self.assertTrue(summary["gb_component_identity_zero"])
        target = self.probe["same_basis_divergence_target"]
        self.assertLess(target["direct_reconstruction_residual_l2"], 1e-12)
        self.assertLess(target["least_squares_residual_l2"], 1e-12)

    def test_gb_lift_rank_test_stays_negative(self) -> None:
        rank_test = self.probe["gb_lift_rank_test"]
        self.assertEqual(rank_test["numeric_rank"], 2)
        self.assertEqual(rank_test["numeric_nullity"], 2)
        self.assertFalse(rank_test["gb_dependence_lifted"])
        self.assertLess(rank_test["null_vector_residual_l2"], 1e-8)

    def test_gatekeepers_and_fingerprint(self) -> None:
        checks = self.payload["gatekeeper_checks"]
        for key in [
            "component_entry_count_12",
            "all_tracefree_sums_zero",
            "gb_component_identity_zero",
            "same_basis_target_reconstructs",
            "gb_lift_rank_is_2",
            "gb_lift_nullity_is_2",
            "full_global_renormalization_not_claimed",
            "no_qw2191_discharge_claimed",
            "no_g1_g3_update_claimed",
            "no_toe_closure_claimed",
        ]:
            self.assertTrue(checks[key], key)
        theorem = self.probe["theorem_export"]
        self.assertEqual(self.probe["theorem_fingerprint_sha256"], sha256_json(theorem))


if __name__ == "__main__":
    unittest.main()
