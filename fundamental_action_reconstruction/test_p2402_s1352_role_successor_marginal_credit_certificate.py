from __future__ import annotations

import hashlib
import json
import subprocess
import sys
import unittest
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parent
SCRIPT = ROOT / "p2402_s1352_role_successor_marginal_credit_certificate.py"
OUT = ROOT / "generated" / "p2402_s1352_role_successor_marginal_credit_certificate.json"
MD = ROOT / "generated" / "p2402_s1352_role_successor_marginal_credit_certificate.md"
PREREQ = ROOT / "p2401_s1351_role_successor_unlock_order_certificate.py"


def sha256_json(obj: Any) -> str:
    return hashlib.sha256(
        json.dumps(obj, sort_keys=True, ensure_ascii=False, separators=(",", ":")).encode("utf-8")
    ).hexdigest()


class P2402RoleSuccessorMarginalCreditCertificateTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        subprocess.run([sys.executable, str(PREREQ)], check=True, cwd=ROOT.parent)
        subprocess.run([sys.executable, str(SCRIPT)], check=True, cwd=ROOT.parent)
        cls.payload = json.loads(OUT.read_text(encoding="utf-8"))
        cls.cert = cls.payload["role_successor_marginal_credit_certificate"]
        cls.theorem = cls.cert["theorem_export"]

    def test_packet_identity_and_outputs(self) -> None:
        self.assertEqual(self.payload["schema_version"], "p2402_s1352_v1")
        self.assertEqual(self.payload["packet_id"], "P2402")
        self.assertEqual(self.payload["stage_id"], "S1352")
        self.assertEqual(self.payload["result_kind"], "ROLE_SUCCESSOR_MARGINAL_CREDIT_CERTIFICATE")
        self.assertTrue(OUT.exists())
        self.assertTrue(MD.exists())

    def test_exact_marginal_credit(self) -> None:
        credit = self.theorem["atom_marginal_credit"]
        self.assertEqual(credit["alpha_geo_electroweak_role_theorem"]["mean_marginal_claim_count"], "11/6")
        self.assertEqual(credit["beta_tors_strict_role_theorem"]["mean_marginal_claim_count"], "4/3")
        self.assertEqual(credit["beta_power_hierarchy_successor_theorem"]["mean_marginal_claim_count"], "5/6")
        self.assertEqual(credit["alpha_geo_electroweak_role_theorem"]["mean_physical_marginal_claim_count"], "3/2")
        self.assertEqual(credit["beta_tors_strict_role_theorem"]["mean_physical_marginal_claim_count"], "1")
        self.assertEqual(credit["beta_power_hierarchy_successor_theorem"]["mean_physical_marginal_claim_count"], "1/2")

    def test_rankings_and_inherited_limits(self) -> None:
        self.assertEqual(self.theorem["order_count"], 6)
        self.assertEqual(self.theorem["total_claim_credit_by_atom_ranking"][0], "alpha_geo_electroweak_role_theorem")
        self.assertEqual(self.theorem["physical_claim_credit_by_atom_ranking"][0], "alpha_geo_electroweak_role_theorem")
        self.assertEqual(self.theorem["p2401_best_order_inherited"], [[
            "alpha_geo_electroweak_role_theorem",
            "beta_tors_strict_role_theorem",
            "beta_power_hierarchy_successor_theorem",
        ]])
        self.assertTrue(self.theorem["p2400_full_role_mask_still_required"])

    def test_fingerprint_gatekeepers_docs_and_limits(self) -> None:
        self.assertEqual(self.cert["theorem_fingerprint_sha256"], sha256_json(self.theorem))
        self.assertTrue(all(self.payload["gatekeeper_checks"].values()))
        self.assertIn("No marginal-credit score exports a role-successor atom.", self.theorem["not_licensed"])
        self.assertIn("role-successor marginal-credit", MD.read_text(encoding="utf-8"))
        self.assertIn("P2402/S1352 role-successor", (ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md").read_text(encoding="utf-8"))
        self.assertIn("P2402/S1352 marginal-credit", (ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md").read_text(encoding="utf-8"))


if __name__ == "__main__":
    unittest.main()
