from __future__ import annotations

import hashlib
import json
import subprocess
import sys
import unittest
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parent
SCRIPT = ROOT / "p2340_s1290_strict_track_b_two_fixture_topological_gluing_compatibility_witness.py"
OUT = ROOT / "generated" / "p2340_s1290_strict_track_b_two_fixture_topological_gluing_compatibility_witness.json"
MD = ROOT / "generated" / "p2340_s1290_strict_track_b_two_fixture_topological_gluing_compatibility_witness.md"


def sha256_json(obj: Any) -> str:
    return hashlib.sha256(
        json.dumps(obj, sort_keys=True, ensure_ascii=False, separators=(",", ":")).encode("utf-8")
    ).hexdigest()


class P2340TrackBTwoFixtureCompatibilityTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        subprocess.run([sys.executable, str(SCRIPT)], check=True, cwd=ROOT.parent)
        cls.payload = json.loads(OUT.read_text(encoding="utf-8"))
        cls.probe = cls.payload["strict_track_b_two_fixture_topological_gluing_compatibility_witness_probe"]
        cls.witness = cls.probe["track_B_two_fixture_witness"]
        cls.compatibility = cls.witness["compatibility"]

    def test_packet_identity_and_outputs(self) -> None:
        self.assertEqual(self.payload["schema_version"], "p2340_s1290_v1")
        self.assertEqual(self.payload["packet_id"], "P2340")
        self.assertEqual(self.payload["stage_id"], "S1290")
        self.assertTrue(OUT.exists())
        self.assertTrue(MD.exists())

    def test_two_fixture_pairings_and_additivity(self) -> None:
        fixtures = self.witness["fixtures"]
        self.assertEqual(fixtures[0]["fixture"], "S^4")
        self.assertEqual(fixtures[0]["normalized_pairing"], "13152087*log(2)/5000000")
        self.assertEqual(fixtures[1]["fixture"], "CP^2")
        self.assertEqual(fixtures[1]["chi"], "3")
        self.assertEqual(fixtures[1]["normalized_pairing"], "39456261*log(2)/10000000")
        self.assertEqual(self.compatibility["chi_disjoint_union"], "5")
        self.assertEqual(self.compatibility["direct_union_pairing"], "13152087*log(2)/2000000")
        self.assertEqual(self.compatibility["additivity_residual"], "0")
        self.assertEqual(self.compatibility["ratio_residual_2cp2_minus_3s4"], "0")

    def test_gatekeepers_and_fingerprint(self) -> None:
        deps = self.probe["current_export_dependencies"]
        self.assertTrue(deps["p2339_s4_fixture_loaded"])
        self.assertTrue(deps["p2339_s4_pairing_matches"])
        for key, value in self.payload["gatekeeper_checks"].items():
            self.assertTrue(value, key)
        theorem = self.probe["theorem_export"]
        self.assertEqual(self.probe["theorem_fingerprint_sha256"], sha256_json(theorem))
        self.assertIn("connected-sum gluing theorem", theorem["not_licensed"])
        self.assertIn("universal boundary/topological normalization over all manifolds", theorem["not_licensed"])


if __name__ == "__main__":
    unittest.main()
