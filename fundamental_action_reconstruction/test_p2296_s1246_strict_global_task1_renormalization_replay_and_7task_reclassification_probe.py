from __future__ import annotations
import json
import subprocess
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent
G = ROOT / "generated"


class TestP2296S1246StrictGlobalTask1RenormalizationReplayAnd7TaskReclassificationProbe(unittest.TestCase):
    def test_packet(self) -> None:
        subprocess.run(
            [sys.executable, str(ROOT / "p2296_s1246_strict_global_task1_renormalization_replay_and_7task_reclassification_probe.py")],
            check=True,
        )
        data = json.loads(
            (G / "p2296_s1246_strict_global_task1_renormalization_replay_and_7task_reclassification_probe.json").read_text(encoding="utf-8")
        )
        self.assertEqual(data["schema_version"], "p2296_s1246_v1")
        probe = data["strict_global_task1_renormalization_replay_and_7task_reclassification_probe"]
        self.assertTrue(probe["symbolic_replay"]["residual_is_zero"])
        self.assertIn("GLOBAL_BIANCHI_EXTENSION_OPEN", probe["task1_reclassification"]["updated_status"])
        self.assertEqual(len(probe["seven_task_status_after_p2296"]), 7)
        g = data["gatekeeper_checks"]
        self.assertTrue(g["p1972_symbolic_residual_zero_replayed"])
        self.assertTrue(g["all_replay_rows_exact_zero"])
        self.assertTrue(g["coefficient_residual_norms_zero"])
        self.assertTrue(g["strict_kernel_moments_positive"])
        self.assertTrue(g["task1_b1_projected_pass"])
        self.assertTrue(g["global_extension_still_open_by_p1985"])
        self.assertTrue(g["no_cutkosky_closure_claimed"])
        self.assertTrue(g["no_background_independence_closure_claimed"])
        self.assertTrue(g["no_selector_closure_claimed"])
        self.assertTrue(g["no_toe_closure_claimed"])


if __name__ == "__main__":
    unittest.main()
