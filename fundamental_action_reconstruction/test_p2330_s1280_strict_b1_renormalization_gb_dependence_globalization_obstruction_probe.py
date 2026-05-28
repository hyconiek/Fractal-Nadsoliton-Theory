from __future__ import annotations

import hashlib
import json
import subprocess
import sys
import unittest
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parent
SCRIPT = ROOT / "p2330_s1280_strict_b1_renormalization_gb_dependence_globalization_obstruction_probe.py"
OUT = ROOT / "generated" / "p2330_s1280_strict_b1_renormalization_gb_dependence_globalization_obstruction_probe.json"
MD = ROOT / "generated" / "p2330_s1280_strict_b1_renormalization_gb_dependence_globalization_obstruction_probe.md"


def sha256_json(obj: Any) -> str:
    return hashlib.sha256(
        json.dumps(obj, sort_keys=True, ensure_ascii=False, separators=(",", ":")).encode("utf-8")
    ).hexdigest()


class P2330B1RenormalizationGBDependenceTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        subprocess.run([sys.executable, str(SCRIPT)], check=True, cwd=ROOT.parent)
        cls.payload = json.loads(OUT.read_text(encoding="utf-8"))
        cls.probe = cls.payload["strict_b1_renormalization_gb_dependence_globalization_obstruction_probe"]

    def test_packet_identity_and_outputs(self) -> None:
        self.assertEqual(self.payload["schema_version"], "p2330_s1280_v1")
        self.assertEqual(self.payload["packet_id"], "P2330")
        self.assertEqual(self.payload["stage_id"], "S1280")
        self.assertTrue(OUT.exists())
        self.assertTrue(MD.exists())

    def test_gb_dependence_certificate(self) -> None:
        cert = self.probe["gb_dependence_certificate"]
        self.assertTrue(cert["symbolic_relation_zero"])
        self.assertEqual(cert["numeric_rank"], 3)
        self.assertEqual(cert["numeric_nullity"], 1)
        self.assertLess(cert["gram_null_residual_l2"], 1e-10)
        self.assertEqual(cert["p1950_reported_rank"], 3)
        self.assertEqual(cert["p1950_reported_nullity"], 1)

    def test_gatekeepers_preserve_no_false_pass(self) -> None:
        checks = self.payload["gatekeeper_checks"]
        for key in [
            "p1848_profiles_loaded",
            "p1950_loaded",
            "p2096_loaded",
            "symbolic_gb_relation_zero",
            "numeric_rank_is_3",
            "numeric_nullity_is_1",
            "null_vector_residual_small",
            "p1950_rank_defect_preserved",
            "p2096_quotient_scope_only_preserved",
            "full_4channel_renormalization_not_claimed",
            "global_background_renormalization_not_claimed",
            "no_qw2191_discharge_claimed",
            "no_g1_g3_update_claimed",
            "no_toe_closure_claimed",
        ]:
            self.assertTrue(checks[key], key)

    def test_theorem_fingerprint(self) -> None:
        theorem = self.probe["theorem_export"]
        self.assertEqual(self.probe["theorem_fingerprint_sha256"], sha256_json(theorem))
        self.assertIn("Gauss-Bonnet", theorem["theorem_name"])
        self.assertFalse(self.probe["quotient_closure_summary"]["full_4channel_independence_proven"])


if __name__ == "__main__":
    unittest.main()
