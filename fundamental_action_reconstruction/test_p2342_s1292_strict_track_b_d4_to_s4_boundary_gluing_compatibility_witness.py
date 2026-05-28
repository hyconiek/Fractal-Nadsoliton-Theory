from __future__ import annotations

import hashlib
import json
import subprocess
import sys
import unittest
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parent
SCRIPT = ROOT / "p2342_s1292_strict_track_b_d4_to_s4_boundary_gluing_compatibility_witness.py"
OUT = ROOT / "generated" / "p2342_s1292_strict_track_b_d4_to_s4_boundary_gluing_compatibility_witness.json"
MD = ROOT / "generated" / "p2342_s1292_strict_track_b_d4_to_s4_boundary_gluing_compatibility_witness.md"


def sha256_json(obj: Any) -> str:
    return hashlib.sha256(
        json.dumps(obj, sort_keys=True, ensure_ascii=False, separators=(",", ":")).encode("utf-8")
    ).hexdigest()


class P2342TrackBD4ToS4BoundaryGluingTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        subprocess.run([sys.executable, str(SCRIPT)], check=True, cwd=ROOT.parent)
        cls.payload = json.loads(OUT.read_text(encoding="utf-8"))
        cls.probe = cls.payload["strict_track_b_d4_to_s4_boundary_gluing_compatibility_witness_probe"]
        cls.witness = cls.probe["track_B_D4_to_S4_boundary_gluing_witness"]
        cls.gluing = cls.witness["gluing_model"]

    def test_packet_identity_and_outputs(self) -> None:
        self.assertEqual(self.payload["schema_version"], "p2342_s1292_v1")
        self.assertEqual(self.payload["packet_id"], "P2342")
        self.assertEqual(self.payload["stage_id"], "S1292")
        self.assertTrue(OUT.exists())
        self.assertTrue(MD.exists())

    def test_d4_to_s4_gluing_values(self) -> None:
        self.assertEqual(self.witness["ledger_coefficient_b_GB"], "13152087*log(2)/(320000000*pi**2)")
        self.assertEqual(self.witness["per_euler_pairing"], "13152087*log(2)/10000000")
        self.assertEqual(self.gluing["chi_glued_from_pieces"], "2")
        self.assertEqual(self.gluing["chi_residual_against_S4"], "0")
        self.assertEqual(self.gluing["glued_topological_number_from_pieces"], "64*pi**2")
        self.assertEqual(self.gluing["S4_topological_number"], "64*pi**2")
        self.assertEqual(self.gluing["gluing_topological_residual"], "0")
        self.assertEqual(self.gluing["left_D4_pairing"], "13152087*log(2)/10000000")
        self.assertEqual(self.gluing["right_D4_pairing"], "13152087*log(2)/10000000")
        self.assertEqual(self.gluing["glued_pairing_from_D4_pieces"], "13152087*log(2)/5000000")
        self.assertEqual(self.gluing["S4_pairing_direct"], "13152087*log(2)/5000000")
        self.assertEqual(self.gluing["pairing_residual_against_S4"], "0")

    def test_dependencies_gatekeepers_and_fingerprint(self) -> None:
        deps = self.probe["current_export_dependencies"]
        self.assertTrue(deps["p2339_s4_fixture_loaded"])
        self.assertTrue(deps["p2341_d4_boundary_fixture_loaded"])
        self.assertTrue(deps["p2339_s4_pairing_matches"])
        self.assertTrue(deps["p2341_d4_pairing_matches"])
        for key, value in self.payload["gatekeeper_checks"].items():
            self.assertTrue(value, key)
        theorem = self.probe["theorem_export"]
        self.assertEqual(self.probe["theorem_fingerprint_sha256"], sha256_json(theorem))
        self.assertIn("general boundary gluing theorem for arbitrary interfaces", theorem["not_licensed"])
        self.assertIn("selector premise or QW-2191 selector discharge", theorem["not_licensed"])
        self.assertIn("ToE closure", theorem["not_licensed"])


if __name__ == "__main__":
    unittest.main()
