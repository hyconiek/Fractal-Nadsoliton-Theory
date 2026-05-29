from __future__ import annotations

import hashlib
import json
import subprocess
import sys
import unittest
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parent
SCRIPT = ROOT / "p2357_s1307_strict_track_b_flat_convex_component_transgression_integration_bridge.py"
OUT = ROOT / "generated" / "p2357_s1307_strict_track_b_flat_convex_component_transgression_integration_bridge.json"
MD = ROOT / "generated" / "p2357_s1307_strict_track_b_flat_convex_component_transgression_integration_bridge.md"


def sha256_json(obj: Any) -> str:
    return hashlib.sha256(
        json.dumps(obj, sort_keys=True, ensure_ascii=False, separators=(",", ":")).encode("utf-8")
    ).hexdigest()


class P2357FlatConvexComponentTransgressionIntegrationBridgeTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        subprocess.run([sys.executable, str(SCRIPT)], check=True, cwd=ROOT.parent)
        cls.payload = json.loads(OUT.read_text(encoding="utf-8"))
        cls.probe = cls.payload["strict_track_b_flat_convex_component_transgression_integration_bridge_probe"]
        cls.bridge = cls.probe["track_B_flat_convex_component_transgression_integration_bridge"]

    def test_packet_identity_and_outputs(self) -> None:
        self.assertEqual(self.payload["schema_version"], "p2357_s1307_v1")
        self.assertEqual(self.payload["packet_id"], "P2357")
        self.assertEqual(self.payload["stage_id"], "S1307")
        self.assertTrue(OUT.exists())
        self.assertTrue(MD.exists())

    def test_local_to_integrated_symbolic_chain(self) -> None:
        self.assertEqual(self.bridge["local_chern_boundary_polynomial"], "8*(K*k1 + K*k2 + K*k3 + 2*k1*k2*k3)")
        self.assertEqual(self.bridge["flat_density"], "16*k1*k2*k3")
        self.assertEqual(self.bridge["flat_density_residual"], "0")
        self.assertEqual(self.bridge["component_integral_rule"], "32*pi**2*eps")
        self.assertEqual(self.bridge["component_integral_residual"], "0")
        self.assertEqual(self.bridge["finite_component_boundary_formula"], "32*pi**2*n")
        self.assertEqual(self.bridge["finite_component_pairing_formula"], "13152087*n*log(2)/10000000")
        self.assertEqual(len(self.bridge["proof_chain"]), 3)

    def test_fixture_bridge_rows_and_replays(self) -> None:
        rows = {row["fixture_id"]: row for row in self.bridge["fixture_bridge_rows"]}
        self.assertEqual(len(rows), 4)
        self.assertTrue(self.bridge["all_bridge_residuals_zero"])
        for row in rows.values():
            self.assertEqual(row["boundary_bridge_residual"], "0", row["fixture_id"])
            self.assertEqual(row["pairing_bridge_residual"], "0", row["fixture_id"])
        self.assertEqual(rows["convex_one_component_degree_plus_one"]["local_transgression_integral_sum"], "32*pi**2")
        self.assertEqual(rows["p2354_round_shell_degree_plus_one_minus_one"]["local_transgression_integral_sum"], "0")
        self.assertEqual(rows["p2355_nonround_shell_degree_plus_one_minus_one"]["local_transgression_integral_sum"], "0")
        self.assertEqual(rows["three_component_oriented_ledger_stress_vector"]["local_transgression_integral_sum"], "-32*pi**2")
        self.assertEqual(self.bridge["p2345_convex_boundary_residual"], "0")
        self.assertEqual(self.bridge["p2348_flat_integral_replay_residual"], "0")
        self.assertTrue(self.bridge["p2356_fixture_residuals_replayed"])
        self.assertEqual(self.bridge["p2356_coverage_rank_replayed"], 4)
        self.assertIn("O3_arbitrary_boundary_transgression_integration", self.bridge["p2353_minimal_cut_replayed"])
        self.assertTrue(self.bridge["o3_cut_partially_attacked_under_hypotheses"])

    def test_dependencies_gatekeepers_and_fingerprint(self) -> None:
        for key, value in self.probe["current_export_dependencies"].items():
            self.assertTrue(value, key)
        for key, value in self.payload["gatekeeper_checks"].items():
            self.assertTrue(value, key)
        theorem = self.probe["theorem_export"]
        self.assertEqual(self.probe["theorem_fingerprint_sha256"], sha256_json(theorem))
        self.assertIn("arbitrary-boundary theorem", theorem["not_licensed"])
        self.assertIn("general transgression theorem without strict convex-component hypotheses", theorem["not_licensed"])
        self.assertIn("selector premise or QW-2191 selector discharge", theorem["not_licensed"])
        self.assertIn("ToE closure", theorem["not_licensed"])


if __name__ == "__main__":
    unittest.main()
