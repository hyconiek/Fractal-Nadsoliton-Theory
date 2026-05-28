from __future__ import annotations

import hashlib
import json
import subprocess
import sys
import unittest
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parent
SCRIPT = ROOT / "p2334_s1284_strict_eom_gb_topological_quotient_scope_theorem.py"
OUT = ROOT / "generated" / "p2334_s1284_strict_eom_gb_topological_quotient_scope_theorem.json"
MD = ROOT / "generated" / "p2334_s1284_strict_eom_gb_topological_quotient_scope_theorem.md"


def sha256_json(obj: Any) -> str:
    return hashlib.sha256(
        json.dumps(obj, sort_keys=True, ensure_ascii=False, separators=(",", ":")).encode("utf-8")
    ).hexdigest()


class P2334EOMGBTopologicalQuotientScopeTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        subprocess.run([sys.executable, str(SCRIPT)], check=True, cwd=ROOT.parent)
        cls.payload = json.loads(OUT.read_text(encoding="utf-8"))
        cls.probe = cls.payload["strict_eom_gb_topological_quotient_scope_probe"]

    def test_packet_identity_and_outputs(self) -> None:
        self.assertEqual(self.payload["schema_version"], "p2334_s1284_v1")
        self.assertEqual(self.payload["packet_id"], "P2334")
        self.assertEqual(self.payload["stage_id"], "S1284")
        self.assertTrue(OUT.exists())
        self.assertTrue(MD.exists())

    def test_universal_eom_quotient_map(self) -> None:
        eom_map = self.probe["universal_eom_quotient_map"]
        self.assertEqual(eom_map["exact_rank"], 2)
        self.assertEqual(eom_map["exact_nullity"], 2)
        self.assertTrue(eom_map["gb_density_relation_maps_to_zero"])
        self.assertTrue(eom_map["pure_gb_eom_variation_maps_to_zero"])
        self.assertEqual(eom_map["max_possible_eom_rank_for_any_component_atlas"], 2)

    def test_numeric_replay_and_existing_exports(self) -> None:
        replay = self.probe["numeric_replay"]
        self.assertLess(replay["gb_density_relation_residual_l2"], 1e-12)
        self.assertLess(replay["pure_gb_eom_relation_residual_l2"], 1e-12)
        audit = self.probe["existing_export_audit"]
        self.assertEqual(audit["p2333_local_spacetime_rank"], 2)
        self.assertEqual(audit["p2333_local_spacetime_nullity"], 2)
        self.assertFalse(audit["genuine_second_tensor_component_atlas_exported_now"])

    def test_gatekeepers_and_fingerprint(self) -> None:
        for key, value in self.payload["gatekeeper_checks"].items():
            self.assertTrue(value, key)
        theorem = self.probe["formal_scope_theorem"]
        self.assertEqual(self.probe["theorem_fingerprint_sha256"], sha256_json(theorem))
        self.assertIn("action-density", theorem["claim"])
        self.assertIn("QW-2191 selector discharge", theorem["not_licensed"])


if __name__ == "__main__":
    unittest.main()
