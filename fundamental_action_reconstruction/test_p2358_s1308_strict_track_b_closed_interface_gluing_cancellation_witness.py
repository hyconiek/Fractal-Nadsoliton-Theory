from __future__ import annotations

import hashlib
import json
import subprocess
import sys
import unittest
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parent
SCRIPT = ROOT / "p2358_s1308_strict_track_b_closed_interface_gluing_cancellation_witness.py"
OUT = ROOT / "generated" / "p2358_s1308_strict_track_b_closed_interface_gluing_cancellation_witness.json"
MD = ROOT / "generated" / "p2358_s1308_strict_track_b_closed_interface_gluing_cancellation_witness.md"


def sha256_json(obj: Any) -> str:
    return hashlib.sha256(
        json.dumps(obj, sort_keys=True, ensure_ascii=False, separators=(",", ":")).encode("utf-8")
    ).hexdigest()


class P2358ClosedInterfaceGluingCancellationTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        subprocess.run([sys.executable, str(SCRIPT)], check=True, cwd=ROOT.parent)
        cls.payload = json.loads(OUT.read_text(encoding="utf-8"))
        cls.probe = cls.payload["strict_track_b_closed_interface_gluing_cancellation_witness_probe"]
        cls.witness = cls.probe["track_B_closed_interface_gluing_cancellation_witness"]

    def test_packet_identity_and_outputs(self) -> None:
        self.assertEqual(self.payload["schema_version"], "p2358_s1308_v1")
        self.assertEqual(self.payload["packet_id"], "P2358")
        self.assertEqual(self.payload["stage_id"], "S1308")
        self.assertTrue(OUT.exists())
        self.assertTrue(MD.exists())

    def test_symbolic_interface_cancellation(self) -> None:
        self.assertEqual(self.witness["target_boundary_per_degree"], "32*pi**2")
        self.assertEqual(self.witness["target_pairing_per_degree"], "13152087*log(2)/10000000")
        self.assertEqual(self.witness["interface_left_boundary"], "32*pi**2")
        self.assertEqual(self.witness["interface_right_boundary"], "-32*pi**2")
        self.assertEqual(self.witness["symbolic_interface_residual"], "0")
        self.assertEqual(self.witness["symbolic_preglue_degree"], "a + b")
        self.assertEqual(self.witness["symbolic_postglue_degree"], "a + b")
        self.assertEqual(self.witness["symbolic_degree_residual"], "0")
        self.assertEqual(self.witness["symbolic_boundary_residual"], "0")
        self.assertEqual(self.witness["symbolic_pairing_residual"], "0")

    def test_gluing_case_rows_and_replays(self) -> None:
        rows = {row["case_id"]: row for row in self.witness["gluing_case_rows"]}
        self.assertEqual(len(rows), 4)
        self.assertTrue(self.witness["all_case_residuals_zero"])
        for row in rows.values():
            self.assertEqual(row["interface_boundary_sum"], "0", row["case_id"])
            self.assertEqual(row["boundary_gluing_residual"], "0", row["case_id"])
            self.assertEqual(row["pairing_gluing_residual"], "0", row["case_id"])
        self.assertEqual(rows["closed_interface_pair_no_external_boundary"]["postglue_boundary_sum"], "0")
        self.assertEqual(rows["single_external_degree_survives_gluing"]["postglue_boundary_sum"], "32*pi**2")
        self.assertEqual(rows["annular_zero_degree_recomposition"]["postglue_boundary_sum"], "0")
        self.assertEqual(rows["nonzero_stress_external_degree_recomposition"]["postglue_boundary_sum"], "32*pi**2")
        self.assertEqual(self.witness["coverage_matrix_rank"], 4)
        self.assertTrue(self.witness["p2356_fixture_residuals_replayed"])
        self.assertTrue(self.witness["p2357_bridge_residuals_replayed"])
        self.assertIn("O5_regularization_corners_and_gluing", self.witness["p2353_minimal_cut_replayed"])
        self.assertTrue(self.witness["o5_cut_partially_attacked_under_closed_interface_hypotheses"])

    def test_dependencies_gatekeepers_and_fingerprint(self) -> None:
        for key, value in self.probe["current_export_dependencies"].items():
            self.assertTrue(value, key)
        for key, value in self.payload["gatekeeper_checks"].items():
            self.assertTrue(value, key)
        theorem = self.probe["theorem_export"]
        self.assertEqual(self.probe["theorem_fingerprint_sha256"], sha256_json(theorem))
        self.assertIn("arbitrary-boundary theorem", theorem["not_licensed"])
        self.assertIn("general gluing theorem", theorem["not_licensed"])
        self.assertIn("corner contribution theorem", theorem["not_licensed"])
        self.assertIn("selector premise or QW-2191 selector discharge", theorem["not_licensed"])
        self.assertIn("ToE closure", theorem["not_licensed"])


if __name__ == "__main__":
    unittest.main()
