from __future__ import annotations

import hashlib
import json
import subprocess
import sys
import unittest
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parent
SCRIPT = ROOT / "p2359_s1309_strict_track_b_s4_cap_complement_gluing_bulk_boundary_witness.py"
OUT = ROOT / "generated" / "p2359_s1309_strict_track_b_s4_cap_complement_gluing_bulk_boundary_witness.json"
MD = ROOT / "generated" / "p2359_s1309_strict_track_b_s4_cap_complement_gluing_bulk_boundary_witness.md"


def sha256_json(obj: Any) -> str:
    return hashlib.sha256(
        json.dumps(obj, sort_keys=True, ensure_ascii=False, separators=(",", ":")).encode("utf-8")
    ).hexdigest()


class P2359S4CapComplementGluingBulkBoundaryTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        subprocess.run([sys.executable, str(SCRIPT)], check=True, cwd=ROOT.parent)
        cls.payload = json.loads(OUT.read_text(encoding="utf-8"))
        cls.probe = cls.payload["strict_track_b_s4_cap_complement_gluing_bulk_boundary_witness_probe"]
        cls.witness = cls.probe["track_B_s4_cap_complement_gluing_bulk_boundary_witness"]

    def test_packet_identity_and_outputs(self) -> None:
        self.assertEqual(self.payload["schema_version"], "p2359_s1309_v1")
        self.assertEqual(self.payload["packet_id"], "P2359")
        self.assertEqual(self.payload["stage_id"], "S1309")
        self.assertTrue(OUT.exists())
        self.assertTrue(MD.exists())

    def test_symbolic_cap_complement_gluing(self) -> None:
        self.assertEqual(self.witness["ledger_coefficient_b_GB"], "13152087*log(2)/(320000000*pi**2)")
        self.assertEqual(self.witness["closed_s4_target_64pi2"], "64*pi**2")
        self.assertEqual(self.witness["north_bulk"], "16*pi**2*(c - 1)**2*(c + 2)")
        self.assertEqual(self.witness["north_boundary"], "-16*pi**2*c*(c**2 - 3)")
        self.assertEqual(self.witness["north_total"], "32*pi**2")
        self.assertEqual(self.witness["south_boundary"], "16*pi**2*c*(c**2 - 3)")
        self.assertEqual(self.witness["south_total"], "32*pi**2")
        self.assertEqual(self.witness["interface_cancellation"], "0")
        self.assertEqual(self.witness["preglue_bulk_boundary_total"], "64*pi**2")
        self.assertEqual(self.witness["postglue_closed_bulk_total"], "64*pi**2")
        self.assertEqual(self.witness["closed_total_residual"], "0")
        self.assertEqual(self.witness["postglue_residual"], "0")
        self.assertEqual(self.witness["gluing_consistency_residual"], "0")
        self.assertEqual(self.witness["c_derivative_residual"], "0")
        self.assertEqual(self.witness["pairing_residual"], "0")

    def test_samples_and_dependency_replays(self) -> None:
        self.assertEqual(len(self.witness["sample_values"]), 3)
        for row in self.witness["sample_values"]:
            self.assertEqual(row["north_total"], "32*pi**2", row["c"])
            self.assertEqual(row["south_total"], "32*pi**2", row["c"])
            self.assertEqual(row["glued_closed_total"], "64*pi**2", row["c"])
        self.assertEqual(self.witness["p2346_bulk_residual"], "0")
        self.assertEqual(self.witness["p2346_boundary_residual"], "0")
        self.assertEqual(self.witness["p2347_boundary_residual"], "0")
        self.assertEqual(self.witness["p2358_symbolic_interface_residual_replayed"], "0")
        self.assertTrue(self.witness["p2358_all_case_residuals_replayed"])
        self.assertIn("O5_regularization_corners_and_gluing", self.witness["p2353_minimal_cut_replayed"])
        self.assertTrue(self.witness["o5_closed_interface_nonflat_bulk_boundary_partially_attacked"])

    def test_dependencies_gatekeepers_and_fingerprint(self) -> None:
        for key, value in self.probe["current_export_dependencies"].items():
            self.assertTrue(value, key)
        for key, value in self.payload["gatekeeper_checks"].items():
            self.assertTrue(value, key)
        theorem = self.probe["theorem_export"]
        self.assertEqual(self.probe["theorem_fingerprint_sha256"], sha256_json(theorem))
        self.assertIn("general Chern-Gauss-Bonnet theorem over arbitrary compact four-manifolds", theorem["not_licensed"])
        self.assertIn("arbitrary-boundary theorem", theorem["not_licensed"])
        self.assertIn("selector premise or QW-2191 selector discharge", theorem["not_licensed"])
        self.assertIn("ToE closure", theorem["not_licensed"])


if __name__ == "__main__":
    unittest.main()
