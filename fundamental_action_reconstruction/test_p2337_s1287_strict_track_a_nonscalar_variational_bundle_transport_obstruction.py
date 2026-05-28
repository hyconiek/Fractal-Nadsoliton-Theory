from __future__ import annotations

import hashlib
import json
import subprocess
import sys
import unittest
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parent
SCRIPT = ROOT / "p2337_s1287_strict_track_a_nonscalar_variational_bundle_transport_obstruction.py"
OUT = ROOT / "generated" / "p2337_s1287_strict_track_a_nonscalar_variational_bundle_transport_obstruction.json"
MD = ROOT / "generated" / "p2337_s1287_strict_track_a_nonscalar_variational_bundle_transport_obstruction.md"


def sha256_json(obj: Any) -> str:
    return hashlib.sha256(
        json.dumps(obj, sort_keys=True, ensure_ascii=False, separators=(",", ":")).encode("utf-8")
    ).hexdigest()


class P2337TrackANonscalarBundleObstructionTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        subprocess.run([sys.executable, str(SCRIPT)], check=True, cwd=ROOT.parent)
        cls.payload = json.loads(OUT.read_text(encoding="utf-8"))
        cls.probe = cls.payload["strict_track_a_nonscalar_variational_bundle_transport_obstruction_probe"]
        cls.obstruction = cls.probe["scalar_lift_obstruction"]

    def test_packet_identity_and_outputs(self) -> None:
        self.assertEqual(self.payload["schema_version"], "p2337_s1287_v1")
        self.assertEqual(self.payload["packet_id"], "P2337")
        self.assertEqual(self.payload["stage_id"], "S1287")
        self.assertTrue(OUT.exists())
        self.assertTrue(MD.exists())

    def test_scalar_lift_obstruction_is_nonzero(self) -> None:
        self.assertEqual(self.probe["component_basis"], ["E_lapse_N", "E_spatial_1", "E_spatial_2", "E_spatial_3"])
        self.assertTrue(self.obstruction["frw_restriction_zero"])
        self.assertTrue(self.obstruction["residual_symbolically_nonzero"])
        self.assertEqual(self.obstruction["residual_term_counts"], [35, 20, 20, 20])
        self.assertGreater(self.obstruction["sample_l2_norm"], 1e-6)
        self.assertTrue(self.obstruction["tracefree_spatial_sum_zero"])

    def test_dependencies_gatekeepers_and_fingerprint(self) -> None:
        deps = self.probe["current_export_dependencies"]
        self.assertTrue(deps["p2333_component_table_loaded"])
        self.assertTrue(deps["p2335_two_track_ledger_loaded"])
        self.assertTrue(deps["p2336_track_a_scalar_transport_loaded"])
        self.assertTrue(deps["p2336_full_tensor_bundle_not_claimed"])
        for key, value in self.payload["gatekeeper_checks"].items():
            self.assertTrue(value, key)
        theorem = self.probe["theorem_export"]
        self.assertEqual(self.probe["theorem_fingerprint_sha256"], sha256_json(theorem))
        self.assertIn("full non-scalar variational-bundle transport theorem", theorem["not_licensed"])
        self.assertIn("transport or normalization of the Track B GB topological ledger", theorem["not_licensed"])


if __name__ == "__main__":
    unittest.main()
