from __future__ import annotations

import hashlib
import json
import subprocess
import sys
import unittest
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parent
SCRIPT = ROOT / "p2345_s1295_strict_track_b_convex_gauss_map_boundary_functional_theorem_and_closure_promise_map.py"
OUT = ROOT / "generated" / "p2345_s1295_strict_track_b_convex_gauss_map_boundary_functional_theorem_and_closure_promise_map.json"
MD = ROOT / "generated" / "p2345_s1295_strict_track_b_convex_gauss_map_boundary_functional_theorem_and_closure_promise_map.md"


def sha256_json(obj: Any) -> str:
    return hashlib.sha256(
        json.dumps(obj, sort_keys=True, ensure_ascii=False, separators=(",", ":")).encode("utf-8")
    ).hexdigest()


class P2345ConvexGaussMapClosurePromiseTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        subprocess.run([sys.executable, str(SCRIPT)], check=True, cwd=ROOT.parent)
        cls.payload = json.loads(OUT.read_text(encoding="utf-8"))
        cls.probe = cls.payload[
            "strict_track_b_convex_gauss_map_boundary_functional_theorem_and_closure_promise_probe"
        ]
        cls.lemma = cls.probe["track_B_convex_gauss_map_boundary_functional_theorem"]
        cls.promise = cls.probe["closure_promise_map"]

    def test_packet_identity_and_outputs(self) -> None:
        self.assertEqual(self.payload["schema_version"], "p2345_s1295_v1")
        self.assertEqual(self.payload["packet_id"], "P2345")
        self.assertEqual(self.payload["stage_id"], "S1295")
        self.assertTrue(OUT.exists())
        self.assertTrue(MD.exists())

    def test_convex_gauss_map_boundary_values(self) -> None:
        self.assertEqual(self.lemma["gauss_map_degree"], "1")
        self.assertEqual(self.lemma["unit_s3_area"], "2*pi**2")
        self.assertEqual(self.lemma["sigma3_integral_convex_class"], "2*pi**2")
        self.assertEqual(self.lemma["boundary_density_normalization_factor"], "16")
        self.assertEqual(self.lemma["boundary_functional_convex_class"], "32*pi**2")
        self.assertEqual(self.lemma["target_topological_number_32pi2_chi"], "32*pi**2")
        self.assertEqual(self.lemma["boundary_residual"], "0")
        self.assertEqual(self.lemma["normalized_track_B_pairing"], "13152087*log(2)/10000000")
        self.assertEqual(self.lemma["normalized_pairing_residual"], "0")
        self.assertEqual(self.lemma["round_P2343_residual_against_lemma"], "0")
        self.assertEqual(self.lemma["spheroidal_P2344_residual_against_lemma"], "0")

    def test_closure_promise_and_reality_readout(self) -> None:
        ranked = self.promise["ranked_best_to_worst"]
        self.assertEqual(ranked[0], "track_B_convex_boundary_functional")
        self.assertGreater(ranked.index("selector_future_choice_lane"), ranked.index("track_B_convex_boundary_functional"))
        readout = self.probe["nadsoliton_reality_readout"]
        self.assertEqual(readout["status"], "interpretive_readout_below_ToE_closure_no_selector_discharge")
        self.assertIn("actual choice of our unique future", readout["not_claimed"])

    def test_dependencies_gatekeepers_and_fingerprint(self) -> None:
        deps = self.probe["current_export_dependencies"]
        self.assertTrue(deps["p2341_boundary_value_rederived"])
        self.assertTrue(deps["p2343_round_case_absorbed"])
        self.assertTrue(deps["p2344_spheroidal_case_absorbed"])
        for key, value in self.payload["gatekeeper_checks"].items():
            self.assertTrue(value, key)
        theorem = self.probe["theorem_export"]
        self.assertEqual(self.probe["theorem_fingerprint_sha256"], sha256_json(theorem))
        self.assertIn("selector premise or QW-2191 selector discharge", theorem["not_licensed"])
        self.assertIn("choice of the unique physical future", theorem["not_licensed"])
        self.assertIn("ToE closure", theorem["not_licensed"])


if __name__ == "__main__":
    unittest.main()
