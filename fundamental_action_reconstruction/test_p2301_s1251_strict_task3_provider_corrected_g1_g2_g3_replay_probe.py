from __future__ import annotations
import hashlib
import json
import subprocess
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent
G = ROOT / "generated"


def sha256_json(payload: object) -> str:
    blob = json.dumps(payload, sort_keys=True, separators=(",", ":"))
    return hashlib.sha256(blob.encode("utf-8")).hexdigest()


class TestP2301S1251StrictTask3ProviderCorrectedG1G2G3ReplayProbe(unittest.TestCase):
    def test_packet(self) -> None:
        subprocess.run(
            [sys.executable, str(ROOT / "p2301_s1251_strict_task3_provider_corrected_g1_g2_g3_replay_probe.py")],
            check=True,
        )
        data = json.loads(
            (G / "p2301_s1251_strict_task3_provider_corrected_g1_g2_g3_replay_probe.json").read_text(encoding="utf-8")
        )
        self.assertEqual(data["schema_version"], "p2301_s1251_v1")
        probe = data["strict_task3_provider_corrected_g1_g2_g3_replay_probe"]
        statuses = [row["status"] for row in probe["recomputed_gap_rows"]]
        self.assertEqual(statuses, ["OPEN", "CLOSED", "OPEN"])
        self.assertAlmostEqual(probe["closure_score"], 1 / 3)
        replay = probe["provider_corrected_transport_replay"]
        self.assertGreater(replay["original_max_transport_residual_l1"], replay["threshold"])
        self.assertLessEqual(replay["provider_corrected_max_transport_residual_l1"], replay["threshold"])
        note = probe["legacy_kernel_context_note"]
        self.assertFalse(note["used_as_mathematical_input"])
        self.assertEqual(sha256_json(probe["theorem_export"]), probe["theorem_fingerprint_sha256"])
        g = data["gatekeeper_checks"]
        self.assertTrue(g["alpha_geo_strict_source_loaded"])
        self.assertTrue(g["alpha_geo_is_four_ln2_not_legacy_import"])
        self.assertTrue(g["p2300_provider_exact"])
        self.assertTrue(g["p2203_transport_rows_loaded"])
        self.assertTrue(g["g1_recomputed_from_p2281"])
        self.assertTrue(g["g2_recomputed_from_p2300_provider_corrected_rows"])
        self.assertTrue(g["g3_recomputed_from_p2280"])
        self.assertTrue(g["only_g2_closed"])
        self.assertTrue(g["closure_score_partial"])
        self.assertTrue(g["task3_not_closed"])
        self.assertTrue(g["no_qw2191_discharge_claimed"])
        self.assertTrue(g["no_selector_closure_claimed"])
        self.assertTrue(g["no_legacy_kernel_role_transfer"])
        self.assertTrue(g["no_toe_closure_claimed"])


if __name__ == "__main__":
    unittest.main()
