from __future__ import annotations

import hashlib
import json
import subprocess
import sys
import unittest
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parent
SCRIPT = ROOT / "p2387_s1337_bathtub_exact_kkt_branch_certificate.py"
OUT = ROOT / "generated" / "p2387_s1337_bathtub_exact_kkt_branch_certificate.json"
MD = ROOT / "generated" / "p2387_s1337_bathtub_exact_kkt_branch_certificate.md"

PREREQ_SCRIPTS = [
    ROOT / "p2384_s1334_symbolic_bathtub_inequality_proof_packet.py",
    ROOT / "p2385_s1335_exact_z12_support_chamber_theorem.py",
    ROOT / "p2386_s1336_bathtub_lp_dual_certificate.py",
]


def sha256_json(obj: Any) -> str:
    return hashlib.sha256(
        json.dumps(obj, sort_keys=True, ensure_ascii=False, separators=(",", ":")).encode("utf-8")
    ).hexdigest()


class P2387BathtubExactKktBranchCertificateTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        for script in PREREQ_SCRIPTS:
            subprocess.run([sys.executable, str(script)], check=True, cwd=ROOT.parent)
        subprocess.run([sys.executable, str(SCRIPT)], check=True, cwd=ROOT.parent)
        cls.payload = json.loads(OUT.read_text(encoding="utf-8"))
        cls.probe = cls.payload["bathtub_exact_kkt_branch_certificate"]
        cls.cert = cls.probe["exact_kkt_branch_certificate"]

    def test_packet_identity_and_outputs(self) -> None:
        self.assertEqual(self.payload["schema_version"], "p2387_s1337_v1")
        self.assertEqual(self.payload["packet_id"], "P2387")
        self.assertEqual(self.payload["stage_id"], "S1337")
        self.assertEqual(self.payload["result_kind"], "BATHTUB_EXACT_KKT_BRANCH_CERTIFICATE")
        self.assertTrue(OUT.exists())
        self.assertTrue(MD.exists())

    def test_exact_branch_certificate(self) -> None:
        self.assertAlmostEqual(self.cert["cap_M"], 1.6, places=15)
        self.assertAlmostEqual(self.cert["cut_t_1_over_M"], 0.625, places=15)
        margins = self.cert["strict_endpoint_margins_from_monotone_q"]
        self.assertTrue(margins["both_positive"])
        self.assertGreater(margins["q_0_minus_lambda"], 0.0)
        self.assertGreater(margins["lambda_minus_q_1"], 0.0)
        branches = self.cert["closed_branch_kkt_certificate"]
        self.assertEqual([row["branch"] for row in branches], ["left_s_less_than_t", "cut_point_measure_zero", "right_s_greater_than_t"])
        self.assertEqual(branches[0]["kkt_products"], {"mu_times_cap_slack": "0", "dual_gap_times_rho": "0"})
        self.assertEqual(branches[2]["kkt_products"], {"mu_times_cap_slack": "0", "dual_gap_times_rho": "0"})

    def test_closed_value_and_audit_grid(self) -> None:
        self.assertLess(self.cert["closed_value_identity"]["absolute_gap"], 1.0e-12)
        audit = self.cert["computational_audit_not_theorem_core"]
        self.assertEqual(audit["grid_points"], 33)
        self.assertLess(audit["max_q_prime_on_grid"], 0.0)
        self.assertGreater(audit["min_left_q_minus_lambda_on_grid_excluding_cut"], 0.0)
        self.assertGreater(audit["min_right_lambda_minus_q_on_grid_excluding_cut"], 0.0)
        self.assertLess(self.probe["p2386_replay"]["packet_value_gap"], 1.0e-12)

    def test_source_status_fingerprint_and_gatekeepers(self) -> None:
        target = self.cert["source_burden_acceptance_target"]
        self.assertEqual(target["status"], "EXACT_KKT_ACCEPTANCE_TARGET_ONLY_SOURCE_STILL_OPEN")
        self.assertIn("rho=M", target["required_branch_saturation"])
        theorem = self.probe["theorem_export"]
        self.assertEqual(self.probe["theorem_fingerprint_sha256"], sha256_json(theorem))
        self.assertIn("strict source theorem deriving rho* or cap M", theorem["not_licensed"])
        self.assertTrue(all(self.payload["gatekeeper_checks"].values()))


if __name__ == "__main__":
    unittest.main()
