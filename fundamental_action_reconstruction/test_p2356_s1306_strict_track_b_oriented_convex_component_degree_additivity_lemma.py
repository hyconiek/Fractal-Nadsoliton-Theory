from __future__ import annotations

import hashlib
import json
import subprocess
import sys
import unittest
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parent
SCRIPT = ROOT / "p2356_s1306_strict_track_b_oriented_convex_component_degree_additivity_lemma.py"
OUT = ROOT / "generated" / "p2356_s1306_strict_track_b_oriented_convex_component_degree_additivity_lemma.json"
MD = ROOT / "generated" / "p2356_s1306_strict_track_b_oriented_convex_component_degree_additivity_lemma.md"


def sha256_json(obj: Any) -> str:
    return hashlib.sha256(
        json.dumps(obj, sort_keys=True, ensure_ascii=False, separators=(",", ":")).encode("utf-8")
    ).hexdigest()


class P2356OrientedConvexComponentDegreeAdditivityTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        subprocess.run([sys.executable, str(SCRIPT)], check=True, cwd=ROOT.parent)
        cls.payload = json.loads(OUT.read_text(encoding="utf-8"))
        cls.probe = cls.payload["strict_track_b_oriented_convex_component_degree_additivity_lemma_probe"]
        cls.lemma = cls.probe["track_B_oriented_convex_component_degree_additivity_lemma"]

    def test_packet_identity_and_outputs(self) -> None:
        self.assertEqual(self.payload["schema_version"], "p2356_s1306_v1")
        self.assertEqual(self.payload["packet_id"], "P2356")
        self.assertEqual(self.payload["stage_id"], "S1306")
        self.assertTrue(OUT.exists())
        self.assertTrue(MD.exists())

    def test_symbolic_formula_and_fixture_residuals(self) -> None:
        self.assertEqual(self.lemma["target_boundary_per_degree"], "32*pi**2")
        self.assertEqual(self.lemma["target_pairing_per_degree"], "13152087*log(2)/10000000")
        self.assertEqual(self.lemma["symbolic_degree_sum_three_component"], "eps1 + eps2 + eps3")
        self.assertEqual(self.lemma["symbolic_boundary_sum_three_component"], "32*pi**2*(eps1 + eps2 + eps3)")
        self.assertEqual(
            self.lemma["symbolic_pairing_sum_three_component"],
            "13152087*(eps1 + eps2 + eps3)*log(2)/10000000",
        )
        self.assertEqual(len(self.lemma["fixture_replays"]), 4)
        self.assertTrue(self.lemma["all_fixture_residuals_zero"])
        for row in self.lemma["fixture_replays"]:
            self.assertEqual(row["boundary_residual"], "0", row["fixture_id"])
            self.assertEqual(row["pairing_residual"], "0", row["fixture_id"])

    def test_fixture_values_and_coverage(self) -> None:
        fixtures = {row["fixture_id"]: row for row in self.lemma["fixture_replays"]}
        self.assertEqual(fixtures["convex_one_component_degree_plus_one"]["boundary_sum"], "32*pi**2")
        self.assertEqual(fixtures["p2354_round_shell_degree_plus_one_minus_one"]["boundary_sum"], "0")
        self.assertEqual(fixtures["p2355_nonround_shell_degree_plus_one_minus_one"]["boundary_sum"], "0")
        self.assertEqual(fixtures["three_component_oriented_ledger_stress_vector"]["boundary_sum"], "-32*pi**2")
        self.assertEqual(fixtures["three_component_oriented_ledger_stress_vector"]["degree_sum"], -1)
        self.assertEqual(self.lemma["coverage_matrix_rank"], 4)
        self.assertEqual(self.lemma["p2354_shell_boundary_residual_replayed"], "0")
        self.assertEqual(self.lemma["p2355_shell_boundary_residual_replayed"], "0")
        self.assertIn("O4_nonconvex_degree_and_orientation_accounting", self.lemma["p2353_minimal_cut_replayed"])
        self.assertTrue(self.lemma["o4_cut_lemma_level_partially_closed_under_hypotheses"])

    def test_dependencies_gatekeepers_and_fingerprint(self) -> None:
        for key, value in self.probe["current_export_dependencies"].items():
            self.assertTrue(value, key)
        for key, value in self.payload["gatekeeper_checks"].items():
            self.assertTrue(value, key)
        theorem = self.probe["theorem_export"]
        self.assertEqual(self.probe["theorem_fingerprint_sha256"], sha256_json(theorem))
        self.assertIn("arbitrary-boundary theorem", theorem["not_licensed"])
        self.assertIn("general nonconvex boundary theorem without convex-component hypotheses", theorem["not_licensed"])
        self.assertIn("selector premise or QW-2191 selector discharge", theorem["not_licensed"])
        self.assertIn("ToE closure", theorem["not_licensed"])


if __name__ == "__main__":
    unittest.main()
