from __future__ import annotations

import hashlib
import json
import subprocess
import sys
import unittest
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parent
SCRIPT = ROOT / "p2386_s1336_bathtub_lp_dual_certificate.py"
OUT = ROOT / "generated" / "p2386_s1336_bathtub_lp_dual_certificate.json"
MD = ROOT / "generated" / "p2386_s1336_bathtub_lp_dual_certificate.md"

PREREQ_SCRIPTS = [
    ROOT / "p2384_s1334_symbolic_bathtub_inequality_proof_packet.py",
    ROOT / "p2385_s1335_exact_z12_support_chamber_theorem.py",
]


def sha256_json(obj: Any) -> str:
    return hashlib.sha256(
        json.dumps(obj, sort_keys=True, ensure_ascii=False, separators=(",", ":")).encode("utf-8")
    ).hexdigest()


class P2386BathtubLpDualCertificateTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        for script in PREREQ_SCRIPTS:
            subprocess.run([sys.executable, str(script)], check=True, cwd=ROOT.parent)
        subprocess.run([sys.executable, str(SCRIPT)], check=True, cwd=ROOT.parent)
        cls.payload = json.loads(OUT.read_text(encoding="utf-8"))
        cls.probe = cls.payload["bathtub_lp_dual_certificate"]
        cls.cert = cls.probe["lp_dual_certificate"]

    def test_packet_identity_and_outputs(self) -> None:
        self.assertEqual(self.payload["schema_version"], "p2386_s1336_v1")
        self.assertEqual(self.payload["packet_id"], "P2386")
        self.assertEqual(self.payload["stage_id"], "S1336")
        self.assertEqual(self.payload["result_kind"], "BATHTUB_LP_DUAL_CERTIFICATE")
        self.assertTrue(OUT.exists())
        self.assertTrue(MD.exists())

    def test_primal_dual_certificate_values(self) -> None:
        self.assertAlmostEqual(self.cert["cap_M"], 1.6, places=15)
        self.assertAlmostEqual(self.cert["early_cut_t_1_over_M"], 0.625, places=15)
        self.assertLess(self.cert["absolute_primal_dual_gap"], 1.0e-12)
        self.assertLess(self.cert["absolute_primal_weight_gap"], 1.0e-12)
        self.assertIn("lambda+mu(s)>=q(s)", self.cert["linear_program"]["dual"])
        self.assertEqual(self.cert["dual_mu_rule"], "mu(s)=max(q(s)-q(1/M),0)")

    def test_dual_feasibility_monotonicity_and_kkt(self) -> None:
        feasibility = self.cert["sampled_dual_feasibility"]
        complementarity = self.cert["sampled_complementarity"]
        self.assertEqual(feasibility["grid_points"], 65)
        self.assertLess(feasibility["max_negative_dual_gap_violation"], 1.0e-12)
        self.assertLess(feasibility["max_q_prime_on_grid"], 0.0)
        self.assertEqual(complementarity["grid_points"], 42)
        self.assertLess(complementarity["max_complementarity_error"], 1.0e-12)
        self.assertIn("q'(s)<0", self.cert["uniqueness_statement"])

    def test_source_burden_fingerprint_and_gatekeepers(self) -> None:
        burden = self.cert["source_burden_translation_for_M_1_6"]
        self.assertEqual(burden["early_interval_length"], 0.625)
        self.assertEqual(burden["early_half_mass"], 0.8)
        self.assertEqual(burden["barycenter"], 0.3125)
        self.assertEqual(burden["shift_from_uniform_barycenter"], 0.1875)
        self.assertEqual(burden["status"], "SOURCE_STILL_OPEN_NON_STRICT_PREMISE_UNTIL_DERIVED")
        theorem = self.probe["theorem_export"]
        self.assertEqual(self.probe["theorem_fingerprint_sha256"], sha256_json(theorem))
        self.assertIn("strict variational source theorem deriving the cap M or the density rho*", theorem["not_licensed"])
        self.assertTrue(all(self.payload["gatekeeper_checks"].values()))


if __name__ == "__main__":
    unittest.main()
