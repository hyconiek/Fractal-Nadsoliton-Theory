from __future__ import annotations

import hashlib
import json
import subprocess
import sys
import unittest
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parent
SCRIPT = ROOT / "p2336_s1286_strict_track_a_second_atlas_finite_part_transport_theorem.py"
OUT = ROOT / "generated" / "p2336_s1286_strict_track_a_second_atlas_finite_part_transport_theorem.json"
MD = ROOT / "generated" / "p2336_s1286_strict_track_a_second_atlas_finite_part_transport_theorem.md"


def sha256_json(obj: Any) -> str:
    return hashlib.sha256(
        json.dumps(obj, sort_keys=True, ensure_ascii=False, separators=(",", ":")).encode("utf-8")
    ).hexdigest()


class P2336TrackASecondAtlasTransportTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        subprocess.run([sys.executable, str(SCRIPT)], check=True, cwd=ROOT.parent)
        cls.payload = json.loads(OUT.read_text(encoding="utf-8"))
        cls.probe = cls.payload["strict_track_a_second_atlas_transport_probe"]
        cls.transport = cls.probe["second_atlas_transport"]

    def test_packet_identity_and_outputs(self) -> None:
        self.assertEqual(self.payload["schema_version"], "p2336_s1286_v1")
        self.assertEqual(self.payload["packet_id"], "P2336")
        self.assertEqual(self.payload["stage_id"], "S1286")
        self.assertTrue(OUT.exists())
        self.assertTrue(MD.exists())

    def test_track_a_transport_is_symbolically_zero_and_invertible(self) -> None:
        self.assertEqual(self.probe["track_A_input"]["basis"], ["E(R2)", "E(Ric2)"])
        self.assertEqual(self.transport["transport_residual_entries"], ["0", "0"])
        self.assertEqual(self.transport["matrix_residual_entries"], ["0", "0", "0", "0"])
        self.assertEqual(self.transport["endpoint_back_transport_residual_entries"], ["0", "0"])
        self.assertEqual(self.transport["determinant_T_A_lambda_1"], "(nu*sigma2 + 1)**2")
        self.assertTrue(self.transport["symbolic_zero_residuals"])

    def test_numeric_replay_and_dependencies(self) -> None:
        replay = self.probe["numeric_transport_replay"]
        self.assertTrue(replay["all_rows_pass"])
        self.assertEqual(len(replay["rows"]), 4)
        deps = self.probe["current_export_dependencies"]
        self.assertTrue(deps["p1973_local_scalar_transport_pass"])
        self.assertTrue(deps["p1973_not_global_background_independence"])
        self.assertTrue(deps["p2334_eom_only_gb_no_go_loaded"])
        self.assertTrue(deps["p2335_two_track_ledger_loaded"])

    def test_gatekeepers_and_fingerprint(self) -> None:
        for key, value in self.payload["gatekeeper_checks"].items():
            self.assertTrue(value, key)
        theorem = self.probe["theorem_export"]
        self.assertEqual(self.probe["theorem_fingerprint_sha256"], sha256_json(theorem))
        self.assertIn("transport or normalization of the Track B GB topological ledger", theorem["not_licensed"])
        self.assertIn("full tensor-component variational bundle transport", theorem["not_licensed"])


if __name__ == "__main__":
    unittest.main()
