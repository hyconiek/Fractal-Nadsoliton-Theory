from __future__ import annotations

import hashlib
import json
import subprocess
import sys
import unittest
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parent
SCRIPT = ROOT / "p2388_s1338_cap_threshold_root_uniqueness_certificate.py"
OUT = ROOT / "generated" / "p2388_s1338_cap_threshold_root_uniqueness_certificate.json"
MD = ROOT / "generated" / "p2388_s1338_cap_threshold_root_uniqueness_certificate.md"

PREREQ_SCRIPTS = [
    ROOT / "p2384_s1334_symbolic_bathtub_inequality_proof_packet.py",
    ROOT / "p2385_s1335_exact_z12_support_chamber_theorem.py",
    ROOT / "p2386_s1336_bathtub_lp_dual_certificate.py",
    ROOT / "p2387_s1337_bathtub_exact_kkt_branch_certificate.py",
]


def sha256_json(obj: Any) -> str:
    return hashlib.sha256(
        json.dumps(obj, sort_keys=True, ensure_ascii=False, separators=(",", ":")).encode("utf-8")
    ).hexdigest()


class P2388CapThresholdRootUniquenessCertificateTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        for script in PREREQ_SCRIPTS:
            subprocess.run([sys.executable, str(script)], check=True, cwd=ROOT.parent)
        subprocess.run([sys.executable, str(SCRIPT)], check=True, cwd=ROOT.parent)
        cls.payload = json.loads(OUT.read_text(encoding="utf-8"))
        cls.probe = cls.payload["cap_threshold_root_uniqueness_certificate"]
        cls.cert = cls.probe["cap_threshold_root_certificate"]

    def test_packet_identity_and_outputs(self) -> None:
        self.assertEqual(self.payload["schema_version"], "p2388_s1338_v1")
        self.assertEqual(self.payload["packet_id"], "P2388")
        self.assertEqual(self.payload["stage_id"], "S1338")
        self.assertEqual(self.payload["result_kind"], "CAP_THRESHOLD_ROOT_UNIQUENESS_CERTIFICATE")
        self.assertTrue(OUT.exists())
        self.assertTrue(MD.exists())

    def test_bracket_bisection_and_newton(self) -> None:
        bracket = self.cert["initial_sign_bracket_around_p2382_value"]
        self.assertTrue(bracket["sign_change"])
        self.assertLess(bracket["F_lo"], 0.0)
        self.assertGreater(bracket["F_hi"], 0.0)
        bisection = self.cert["bisection_certificate"]
        self.assertLessEqual(bisection["width"], 1.0e-15)
        self.assertAlmostEqual(bisection["mid"], 1.574821357435363, places=12)
        self.assertLess(self.cert["newton_replay"]["final_abs_residual"], 1.0e-14)

    def test_derivative_uniqueness_and_m16_acceptance(self) -> None:
        audit = self.cert["derivative_audit"]
        self.assertTrue(audit["all_sampled_positive"])
        self.assertIn("P2384 cap-derivative sign proof", audit["proof_role"])
        accepted = self.cert["accepted_cap_replay"]
        self.assertTrue(accepted["accepted"])
        self.assertAlmostEqual(accepted["M_test"], 1.6)
        self.assertGreater(accepted["M_test_minus_root_mid"], 0.025)
        self.assertGreater(accepted["F_M_test"], 0.0)

    def test_source_burden_fingerprint_and_gatekeepers(self) -> None:
        burden = self.cert["source_burden_translation"]
        self.assertEqual(burden["status"], "UNIQUE_THRESHOLD_ACCEPTANCE_TARGET_ONLY_SOURCE_STILL_OPEN")
        self.assertAlmostEqual(burden["M_1_6_early_interval_length"], 0.625)
        self.assertAlmostEqual(burden["M_1_6_early_half_mass"], 0.8)
        self.assertAlmostEqual(burden["M_1_6_barycenter"], 0.3125)
        self.assertAlmostEqual(burden["M_1_6_barycenter_shift_from_uniform"], 0.1875)
        theorem = self.probe["theorem_export"]
        self.assertEqual(self.probe["theorem_fingerprint_sha256"], sha256_json(theorem))
        self.assertIn("strict source theorem deriving M or the front-loaded density", theorem["not_licensed"])
        self.assertTrue(all(self.payload["gatekeeper_checks"].values()))


if __name__ == "__main__":
    unittest.main()
