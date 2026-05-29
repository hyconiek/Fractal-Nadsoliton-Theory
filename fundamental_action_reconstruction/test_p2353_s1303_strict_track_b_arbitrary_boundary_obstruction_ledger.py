from __future__ import annotations

import hashlib
import json
import subprocess
import sys
import unittest
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parent
SCRIPT = ROOT / "p2353_s1303_strict_track_b_arbitrary_boundary_obstruction_ledger.py"
OUT = ROOT / "generated" / "p2353_s1303_strict_track_b_arbitrary_boundary_obstruction_ledger.json"
MD = ROOT / "generated" / "p2353_s1303_strict_track_b_arbitrary_boundary_obstruction_ledger.md"


def sha256_json(obj: Any) -> str:
    return hashlib.sha256(
        json.dumps(obj, sort_keys=True, ensure_ascii=False, separators=(",", ":")).encode("utf-8")
    ).hexdigest()


class P2353ArbitraryBoundaryObstructionLedgerTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        subprocess.run([sys.executable, str(SCRIPT)], check=True, cwd=ROOT.parent)
        cls.payload = json.loads(OUT.read_text(encoding="utf-8"))
        cls.probe = cls.payload["strict_track_b_arbitrary_boundary_obstruction_ledger_probe"]
        cls.ledger = cls.probe["track_B_arbitrary_boundary_obstruction_ledger"]

    def test_packet_identity_and_outputs(self) -> None:
        self.assertEqual(self.payload["schema_version"], "p2353_s1303_v1")
        self.assertEqual(self.payload["packet_id"], "P2353")
        self.assertEqual(self.payload["stage_id"], "S1303")
        self.assertTrue(OUT.exists())
        self.assertTrue(MD.exists())

    def test_replay_residuals_and_target_pairing(self) -> None:
        self.assertEqual(self.ledger["ledger_coefficient_b_GB"], "13152087*log(2)/(320000000*pi**2)")
        self.assertEqual(self.ledger["target_boundary_functional"], "32*pi**2")
        self.assertEqual(self.ledger["target_pairing_per_euler_number"], "13152087*log(2)/10000000")
        residuals = self.ledger["replay_residuals"]
        self.assertEqual(residuals["p2345_convex_gauss_map_boundary_minus_32pi2"], "0")
        self.assertEqual(residuals["p2348_flat_chern_polynomial_replay_minus_32pi2"], "0")
        self.assertEqual(residuals["p2352_boundary_residual_vector"], ["0", "0", "0", "0", "0"])
        self.assertEqual(residuals["p2352_pairing_residual_vector"], ["0", "0", "0", "0", "0"])

    def test_obstruction_counts_vectors_and_cut(self) -> None:
        self.assertEqual(self.ledger["required_obligation_count"], 6)
        self.assertEqual(self.ledger["discharged_required_count"], 2)
        self.assertEqual(self.ledger["partial_required_count"], 1)
        self.assertEqual(self.ledger["open_required_count"], 4)
        self.assertEqual(self.ledger["hard_open_required_count"], 3)
        self.assertEqual(self.ledger["closure_score_discharged_only"], "1/3")
        self.assertEqual(self.ledger["closure_score_with_partial_half_credit"], "5/12")
        self.assertEqual(self.ledger["required_obstruction_vector"], [0, 0, 1, 1, 1, 1])
        self.assertEqual(self.ledger["hard_required_obstruction_vector"], [0, 0, 1, 1, 1, 0])
        self.assertEqual(
            self.ledger["minimal_next_missing_cut"],
            [
                "O3_arbitrary_boundary_transgression_integration",
                "O4_nonconvex_degree_and_orientation_accounting",
                "O5_regularization_corners_and_gluing",
            ],
        )
        self.assertGreaterEqual(self.ledger["support_matrix_rank"], 4)

    def test_dependencies_gatekeepers_and_fingerprint(self) -> None:
        for key, value in self.probe["current_export_dependencies"].items():
            self.assertTrue(value, key)
        for key, value in self.payload["gatekeeper_checks"].items():
            self.assertTrue(value, key)
        theorem = self.probe["theorem_export"]
        self.assertEqual(self.probe["theorem_fingerprint_sha256"], sha256_json(theorem))
        self.assertIn("arbitrary-boundary theorem", theorem["not_licensed"])
        self.assertIn("selector premise or QW-2191 selector discharge", theorem["not_licensed"])
        self.assertIn("ToE closure", theorem["not_licensed"])
        self.assertIn("Do not present another controlled sample as universal closure", theorem["next_honest_step"])


if __name__ == "__main__":
    unittest.main()
