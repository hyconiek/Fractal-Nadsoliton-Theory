from __future__ import annotations

import hashlib
import json
import subprocess
import sys
import unittest
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parent
SCRIPT = ROOT / "p2335_s1285_strict_renormalization_two_track_non_gb_eom_and_gb_topological_ledger.py"
OUT = ROOT / "generated" / "p2335_s1285_strict_renormalization_two_track_non_gb_eom_and_gb_topological_ledger.json"
MD = ROOT / "generated" / "p2335_s1285_strict_renormalization_two_track_non_gb_eom_and_gb_topological_ledger.md"


def sha256_json(obj: Any) -> str:
    return hashlib.sha256(
        json.dumps(obj, sort_keys=True, ensure_ascii=False, separators=(",", ":")).encode("utf-8")
    ).hexdigest()


class P2335TwoTrackRenormalizationLedgerTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        subprocess.run([sys.executable, str(SCRIPT)], check=True, cwd=ROOT.parent)
        cls.payload = json.loads(OUT.read_text(encoding="utf-8"))
        cls.probe = cls.payload["strict_two_track_renormalization_ledger_probe"]
        cls.ledger = cls.probe["two_track_ledger"]

    def test_packet_identity_and_outputs(self) -> None:
        self.assertEqual(self.payload["schema_version"], "p2335_s1285_v1")
        self.assertEqual(self.payload["packet_id"], "P2335")
        self.assertEqual(self.payload["stage_id"], "S1285")
        self.assertTrue(OUT.exists())
        self.assertTrue(MD.exists())

    def test_two_track_split_reconstructs_action_density(self) -> None:
        split = self.ledger["exact_action_density_split"]
        self.assertTrue(split["action_density_identity_zero"])
        self.assertEqual(split["reconstruction_residuals"], {"R2": "0", "Ric2": "0", "Riem2": "0"})
        self.assertEqual(split["b_R2_equals"], "a_R2 - a_Riem2")
        self.assertEqual(split["b_GB_equals"], "a_GB + a_Riem2")

    def test_tracks_are_separated(self) -> None:
        track_a = self.ledger["track_A_non_gb_eom_transportable_quotient"]
        track_b = self.ledger["track_B_gb_topological_counterterm_ledger"]
        self.assertEqual(track_a["basis"], ["E(R2)", "E(Ric2)"])
        self.assertEqual(len(track_a["numeric_coefficients"]), 2)
        self.assertEqual(track_b["status"], "separate_topological_counterterm_ledger_not_EOM_transportable")
        self.assertIn("boundary/topological-number", track_b["requires_future_witness"])

    def test_dependencies_gatekeepers_and_fingerprint(self) -> None:
        deps = self.probe["current_export_dependencies"]
        self.assertTrue(deps["p2334_current_exports_quotient_scope_only"])
        for key, value in self.payload["gatekeeper_checks"].items():
            self.assertTrue(value, key)
        theorem = self.probe["theorem_export"]
        self.assertEqual(self.probe["theorem_fingerprint_sha256"], sha256_json(theorem))
        self.assertIn("independent a_GB measurement", theorem["not_licensed"])


if __name__ == "__main__":
    unittest.main()
