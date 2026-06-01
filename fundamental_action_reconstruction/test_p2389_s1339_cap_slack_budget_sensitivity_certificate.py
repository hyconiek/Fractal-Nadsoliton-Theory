from __future__ import annotations

import hashlib
import json
import subprocess
import sys
import unittest
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parent
SCRIPT = ROOT / "p2389_s1339_cap_slack_budget_sensitivity_certificate.py"
OUT = ROOT / "generated" / "p2389_s1339_cap_slack_budget_sensitivity_certificate.json"
MD = ROOT / "generated" / "p2389_s1339_cap_slack_budget_sensitivity_certificate.md"

PREREQ_SCRIPTS = [
    ROOT / "p2384_s1334_symbolic_bathtub_inequality_proof_packet.py",
    ROOT / "p2385_s1335_exact_z12_support_chamber_theorem.py",
    ROOT / "p2386_s1336_bathtub_lp_dual_certificate.py",
    ROOT / "p2387_s1337_bathtub_exact_kkt_branch_certificate.py",
    ROOT / "p2388_s1338_cap_threshold_root_uniqueness_certificate.py",
]


def sha256_json(obj: Any) -> str:
    return hashlib.sha256(
        json.dumps(obj, sort_keys=True, ensure_ascii=False, separators=(",", ":")).encode("utf-8")
    ).hexdigest()


class P2389CapSlackBudgetSensitivityCertificateTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        for script in PREREQ_SCRIPTS:
            subprocess.run([sys.executable, str(script)], check=True, cwd=ROOT.parent)
        subprocess.run([sys.executable, str(SCRIPT)], check=True, cwd=ROOT.parent)
        cls.payload = json.loads(OUT.read_text(encoding="utf-8"))
        cls.probe = cls.payload["cap_slack_budget_sensitivity_certificate"]
        cls.cert = cls.probe["cap_slack_budget_certificate"]
        cls.budget = cls.cert["derivative_and_margin_budget"]
        cls.geometry = cls.cert["source_geometry_budget"]

    def test_packet_identity_and_outputs(self) -> None:
        self.assertEqual(self.payload["schema_version"], "p2389_s1339_v1")
        self.assertEqual(self.payload["packet_id"], "P2389")
        self.assertEqual(self.payload["stage_id"], "S1339")
        self.assertEqual(self.payload["result_kind"], "CAP_SLACK_BUDGET_SENSITIVITY_CERTIFICATE")
        self.assertTrue(OUT.exists())
        self.assertTrue(MD.exists())

    def test_margin_derivative_and_sensitivity_budget(self) -> None:
        self.assertGreater(self.budget["cap_surplus_M_minus_root"], 0.025)
        self.assertAlmostEqual(self.budget["F_M_margin"], 0.02569532538409236)
        self.assertGreater(self.budget["min_derivative_grid"]["F_prime"], 0.0)
        self.assertTrue(self.budget["mean_value_slope_inside_grid_bounds"])
        sensitivity = self.budget["root_sensitivity_to_additive_threshold_error_interval"]
        self.assertGreater(sensitivity["dM_dT_upper_1_over_min_F_prime"], sensitivity["dM_dT_lower_1_over_max_F_prime"])
        self.assertGreater(self.budget["full_margin_can_absorb_additive_target_error_below"], 0.025)

    def test_source_geometry_slack_budget(self) -> None:
        self.assertEqual(self.geometry["status"], "SLACK_BUDGET_ONLY_SOURCE_STILL_OPEN")
        self.assertAlmostEqual(self.geometry["accepted_cap_M"], 1.6)
        self.assertGreater(self.geometry["early_interval_shortening_vs_threshold"], 0.0)
        self.assertGreater(self.geometry["early_half_mass_surplus_vs_threshold"], 0.0)
        self.assertGreater(self.geometry["barycenter_left_shift_vs_threshold"], 0.0)
        self.assertGreater(self.geometry["uniform_shift_surplus_vs_threshold"], 0.0)

    def test_replay_fingerprint_and_gatekeepers(self) -> None:
        self.assertAlmostEqual(self.probe["p2388_replay"]["p2388_F_M_test"], self.budget["F_M_margin"])
        theorem = self.probe["theorem_export"]
        self.assertEqual(self.probe["theorem_fingerprint_sha256"], sha256_json(theorem))
        self.assertIn("strict source theorem deriving M=1.6 or the front-loaded density", theorem["not_licensed"])
        self.assertTrue(all(self.payload["gatekeeper_checks"].values()))


if __name__ == "__main__":
    unittest.main()
