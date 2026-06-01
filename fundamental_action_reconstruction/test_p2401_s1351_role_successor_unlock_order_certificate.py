from __future__ import annotations

import hashlib
import json
import subprocess
import sys
import unittest
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parent
SCRIPT = ROOT / "p2401_s1351_role_successor_unlock_order_certificate.py"
OUT = ROOT / "generated" / "p2401_s1351_role_successor_unlock_order_certificate.json"
MD = ROOT / "generated" / "p2401_s1351_role_successor_unlock_order_certificate.md"
PREREQ = ROOT / "p2400_s1350_nearest_lift_role_successor_lattice_certificate.py"


def sha256_json(obj: Any) -> str:
    return hashlib.sha256(
        json.dumps(obj, sort_keys=True, ensure_ascii=False, separators=(",", ":")).encode("utf-8")
    ).hexdigest()


class P2401RoleSuccessorUnlockOrderCertificateTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        subprocess.run([sys.executable, str(PREREQ)], check=True, cwd=ROOT.parent)
        subprocess.run([sys.executable, str(SCRIPT)], check=True, cwd=ROOT.parent)
        cls.payload = json.loads(OUT.read_text(encoding="utf-8"))
        cls.cert = cls.payload["role_successor_unlock_order_certificate"]
        cls.theorem = cls.cert["theorem_export"]

    def test_packet_identity_and_outputs(self) -> None:
        self.assertEqual(self.payload["schema_version"], "p2401_s1351_v1")
        self.assertEqual(self.payload["packet_id"], "P2401")
        self.assertEqual(self.payload["stage_id"], "S1351")
        self.assertEqual(self.payload["result_kind"], "ROLE_SUCCESSOR_UNLOCK_ORDER_CERTIFICATE")
        self.assertTrue(OUT.exists())
        self.assertTrue(MD.exists())

    def test_best_order_and_expected_steps(self) -> None:
        self.assertEqual(self.theorem["permutation_count"], 6)
        self.assertEqual(self.theorem["best_total_unlock_step_sum"], 9)
        self.assertEqual(self.theorem["best_orders_by_total_unlock_step_sum"], [[
            "alpha_geo_electroweak_role_theorem",
            "beta_tors_strict_role_theorem",
            "beta_power_hierarchy_successor_theorem",
        ]])
        self.assertEqual(self.theorem["expected_unlock_step_by_claim_over_uniform_orders"]["legacy_weinberg_sin2_theta_role_transfer"], "2")
        self.assertEqual(self.theorem["expected_unlock_step_by_claim_over_uniform_orders"]["legacy_alpha_em_inverse_role_transfer"], "8/3")

    def test_prefix_unlock_limits(self) -> None:
        summary = self.theorem["first_atom_summary"]
        self.assertEqual(summary["alpha_geo_electroweak_role_theorem"]["step1_unlocked_claims_union"], ["legacy_weinberg_sin2_theta_role_transfer"])
        self.assertEqual(summary["beta_tors_strict_role_theorem"]["step1_unlocked_claims_union"], [])
        self.assertEqual(summary["beta_power_hierarchy_successor_theorem"]["step1_unlocked_claims_union"], [])
        rows = self.cert["unlock_order_enumeration"]["rows"]
        self.assertTrue(all(not row["prefixes"][2]["all_roles_and_toe_conditionally_ready"] for row in rows))
        self.assertTrue(all(row["prefixes"][3]["all_roles_and_toe_conditionally_ready"] for row in rows))

    def test_fingerprint_gatekeepers_docs_and_limits(self) -> None:
        self.assertEqual(self.cert["theorem_fingerprint_sha256"], sha256_json(self.theorem))
        self.assertTrue(all(self.payload["gatekeeper_checks"].values()))
        self.assertTrue(self.theorem["p2400_only_full_role_mask_closes_toe"])
        self.assertIn("No one-atom or two-atom prefix licenses role transfer or ToE closure.", self.theorem["not_licensed"])
        self.assertIn("role-successor unlock-order", MD.read_text(encoding="utf-8"))
        self.assertIn("P2401/S1351 role-successor", (ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md").read_text(encoding="utf-8"))
        self.assertIn("P2401/S1351 unlock-order", (ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md").read_text(encoding="utf-8"))


if __name__ == "__main__":
    unittest.main()
