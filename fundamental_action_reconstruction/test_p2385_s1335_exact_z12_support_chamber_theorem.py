from __future__ import annotations

import hashlib
import json
import subprocess
import sys
import unittest
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parent
SCRIPT = ROOT / "p2385_s1335_exact_z12_support_chamber_theorem.py"
OUT = ROOT / "generated" / "p2385_s1335_exact_z12_support_chamber_theorem.json"
MD = ROOT / "generated" / "p2385_s1335_exact_z12_support_chamber_theorem.md"

PREREQ_SCRIPTS = [
    ROOT / "p2384_s1334_symbolic_bathtub_inequality_proof_packet.py",
]

EXPECTED_PAIR_DISTRIBUTION = {
    "0,0": 12,
    "0,2": 12,
    "0,4": 12,
    "1,1": 72,
    "1,2": 96,
    "1,3": 72,
    "2,0": 12,
    "2,1": 96,
    "2,2": 204,
    "2,3": 48,
    "3,1": 72,
    "3,2": 48,
    "3,3": 24,
    "4,0": 12,
}


def sha256_json(obj: Any) -> str:
    return hashlib.sha256(
        json.dumps(obj, sort_keys=True, ensure_ascii=False, separators=(",", ":")).encode("utf-8")
    ).hexdigest()


class P2385ExactZ12SupportChamberTheoremTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        for script in PREREQ_SCRIPTS:
            subprocess.run([sys.executable, str(script)], check=True, cwd=ROOT.parent)
        subprocess.run([sys.executable, str(SCRIPT)], check=True, cwd=ROOT.parent)
        cls.payload = json.loads(OUT.read_text(encoding="utf-8"))
        cls.probe = cls.payload["exact_z12_support_chamber_theorem"]
        cls.chamber = cls.probe["exact_support_chamber_certificate"]
        cls.replay = cls.probe["p2382_cap_chamber_replay_summary"]

    def test_packet_identity_and_outputs(self) -> None:
        self.assertEqual(self.payload["schema_version"], "p2385_s1335_v1")
        self.assertEqual(self.payload["packet_id"], "P2385")
        self.assertEqual(self.payload["stage_id"], "S1335")
        self.assertEqual(self.payload["result_kind"], "EXACT_Z12_SUPPORT_CHAMBER_THEOREM")
        self.assertTrue(OUT.exists())
        self.assertTrue(MD.exists())

    def test_pair_distribution_and_d5_path_targets(self) -> None:
        self.assertEqual(self.chamber["total_supports_checked"], 792)
        self.assertEqual(self.chamber["pair_distribution"], EXPECTED_PAIR_DISTRIBUTION)
        self.assertEqual(self.chamber["max_h5"], 4)
        self.assertEqual(self.chamber["max_h5_pair_distribution"], {"0,4": 12})
        self.assertEqual(self.chamber["target_count"], 12)
        self.assertEqual(self.chamber["target_pair_distribution"], {"0,4": 12})
        self.assertTrue(self.chamber["d5_path_supports_match_enumerated_targets"])

    def test_integer_gap_chamber_and_boundary_tie(self) -> None:
        self.assertEqual(self.chamber["invalid_gap_rows"], [])
        self.assertEqual(len(self.chamber["boundary_ties_at_r_1_3"]), 1)
        tie = self.chamber["boundary_ties_at_r_1_3"][0]
        self.assertEqual((tie["h1"], tie["h5"], tie["count"]), (3, 3, 24))
        non_target_rows = [row for row in self.chamber["gap_table"] if not row["target_pair_0_4"]]
        self.assertTrue(all(row["score_gap_from_0_4_over_b_at_r_equals_1_3_numerator"] >= 0 for row in non_target_rows))
        self.assertIn("r<1/3", self.chamber["proof_rule"])

    def test_p2382_replay_fingerprint_and_gatekeepers(self) -> None:
        self.assertTrue(self.replay["cap_test_d5_chamber"])
        self.assertEqual(self.replay["cap_test_pair_distribution"], {"0,4": 12})
        self.assertFalse(self.replay["below_threshold_d5_chamber"])
        theorem = self.probe["theorem_export"]
        self.assertEqual(self.probe["theorem_fingerprint_sha256"], sha256_json(theorem))
        self.assertIn("exact finite support theorem", theorem["claim"])
        self.assertIn("strict source theorem deriving the cap M or bang-bang density", theorem["not_licensed"])
        self.assertTrue(all(self.payload["gatekeeper_checks"].values()))


if __name__ == "__main__":
    unittest.main()
