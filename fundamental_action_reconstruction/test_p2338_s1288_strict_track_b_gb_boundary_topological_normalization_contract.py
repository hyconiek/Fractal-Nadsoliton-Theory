from __future__ import annotations

import hashlib
import json
import subprocess
import sys
import unittest
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parent
SCRIPT = ROOT / "p2338_s1288_strict_track_b_gb_boundary_topological_normalization_contract.py"
OUT = ROOT / "generated" / "p2338_s1288_strict_track_b_gb_boundary_topological_normalization_contract.json"
MD = ROOT / "generated" / "p2338_s1288_strict_track_b_gb_boundary_topological_normalization_contract.md"


def sha256_json(obj: Any) -> str:
    return hashlib.sha256(
        json.dumps(obj, sort_keys=True, ensure_ascii=False, separators=(",", ":")).encode("utf-8")
    ).hexdigest()


class P2338TrackBGBBoundaryTopologicalContractTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        subprocess.run([sys.executable, str(SCRIPT)], check=True, cwd=ROOT.parent)
        cls.payload = json.loads(OUT.read_text(encoding="utf-8"))
        cls.probe = cls.payload["strict_track_b_gb_boundary_topological_normalization_contract_probe"]
        cls.contract = cls.probe["track_B_contract"]

    def test_packet_identity_and_outputs(self) -> None:
        self.assertEqual(self.payload["schema_version"], "p2338_s1288_v1")
        self.assertEqual(self.payload["packet_id"], "P2338")
        self.assertEqual(self.payload["stage_id"], "S1288")
        self.assertTrue(OUT.exists())
        self.assertTrue(MD.exists())

    def test_contract_coefficients_and_missing_fields(self) -> None:
        self.assertEqual(self.contract["track"], "Track B GB topological counterterm ledger only")
        self.assertEqual(self.contract["ledger_coefficient_b_GB"], "13152087*log(2)/(320000000*pi**2)")
        self.assertEqual(self.contract["conditional_per_euler_number_pairing"], "13152087*log(2)/10000000")
        self.assertFalse(self.contract["current_exports_supply_required_fields"])
        self.assertEqual(len(self.contract["required_future_witness_fields"]), 5)
        self.assertIn("boundary_correction_functional", self.contract["required_future_witness_fields"])

    def test_dependencies_gatekeepers_and_fingerprint(self) -> None:
        deps = self.probe["current_export_dependencies"]
        self.assertTrue(deps["p2334_eom_gb_no_go_loaded"])
        self.assertTrue(deps["p2335_two_track_ledger_loaded"])
        self.assertTrue(deps["p2337_track_a_obstruction_loaded"])
        for key, value in self.payload["gatekeeper_checks"].items():
            self.assertTrue(value, key)
        theorem = self.probe["theorem_export"]
        self.assertEqual(self.probe["theorem_fingerprint_sha256"], sha256_json(theorem))
        self.assertIn("boundary/topological-number normalization theorem", theorem["not_licensed"])
        self.assertIn("independent a_GB measurement", theorem["not_licensed"])


if __name__ == "__main__":
    unittest.main()
