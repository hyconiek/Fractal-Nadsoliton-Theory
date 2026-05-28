from __future__ import annotations

import hashlib
import json
import subprocess
import sys
import unittest
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parent
SCRIPT = ROOT / "p2339_s1289_strict_track_b_s4_boundaryless_topological_normalization_witness.py"
OUT = ROOT / "generated" / "p2339_s1289_strict_track_b_s4_boundaryless_topological_normalization_witness.json"
MD = ROOT / "generated" / "p2339_s1289_strict_track_b_s4_boundaryless_topological_normalization_witness.md"


def sha256_json(obj: Any) -> str:
    return hashlib.sha256(
        json.dumps(obj, sort_keys=True, ensure_ascii=False, separators=(",", ":")).encode("utf-8")
    ).hexdigest()


class P2339TrackBS4NormalizationWitnessTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        subprocess.run([sys.executable, str(SCRIPT)], check=True, cwd=ROOT.parent)
        cls.payload = json.loads(OUT.read_text(encoding="utf-8"))
        cls.probe = cls.payload["strict_track_b_s4_boundaryless_topological_normalization_witness_probe"]
        cls.witness = cls.probe["track_B_S4_witness"]

    def test_packet_identity_and_outputs(self) -> None:
        self.assertEqual(self.payload["schema_version"], "p2339_s1289_v1")
        self.assertEqual(self.payload["packet_id"], "P2339")
        self.assertEqual(self.payload["stage_id"], "S1289")
        self.assertTrue(OUT.exists())
        self.assertTrue(MD.exists())

    def test_s4_fixture_pairing(self) -> None:
        self.assertTrue(self.witness["supplies_all_p2338_required_fields"])
        self.assertEqual(self.witness["chi_S4"], "2")
        self.assertEqual(self.witness["boundary_correction_value"], "0")
        self.assertEqual(self.witness["gb_topological_number_32pi2_chi"], "64*pi**2")
        self.assertEqual(self.witness["normalized_S4_pairing"], "13152087*log(2)/5000000")
        self.assertEqual(self.witness["pairing_residual"], "0")

    def test_gatekeepers_and_fingerprint(self) -> None:
        deps = self.probe["current_export_dependencies"]
        self.assertTrue(deps["p2335_two_track_ledger_loaded"])
        self.assertTrue(deps["p2338_contract_loaded"])
        for key, value in self.payload["gatekeeper_checks"].items():
            self.assertTrue(value, key)
        theorem = self.probe["theorem_export"]
        self.assertEqual(self.probe["theorem_fingerprint_sha256"], sha256_json(theorem))
        self.assertIn("independent a_GB measurement separate from the P2335 ledger coefficient", theorem["not_licensed"])
        self.assertIn("universal boundary/topological normalization over all manifolds", theorem["not_licensed"])


if __name__ == "__main__":
    unittest.main()
