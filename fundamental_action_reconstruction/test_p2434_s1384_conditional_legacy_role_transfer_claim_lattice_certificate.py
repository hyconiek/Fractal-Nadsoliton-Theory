import json
import subprocess
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent
SCRIPT = ROOT / "p2434_s1384_conditional_legacy_role_transfer_claim_lattice_certificate.py"
OUT = ROOT / "generated" / "p2434_s1384_conditional_legacy_role_transfer_claim_lattice_certificate.json"
MD = ROOT / "generated" / "p2434_s1384_conditional_legacy_role_transfer_claim_lattice_certificate.md"
P2433 = ROOT / "generated" / "p2433_s1383_source_selector_convergence_role_transfer_gate_certificate.json"


class P2434ConditionalLegacyRoleTransferClaimLatticeCertificateTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        if not P2433.exists():
            subprocess.run([sys.executable, str(ROOT / "p2433_s1383_source_selector_convergence_role_transfer_gate_certificate.py")], check=True)
        subprocess.run([sys.executable, str(SCRIPT)], check=True)
        cls.payload = json.loads(OUT.read_text(encoding="utf-8"))
        cls.cert = cls.payload["conditional_legacy_role_transfer_claim_lattice_certificate"]
        cls.theorem = cls.cert["theorem_export"]

    def test_identity_outputs_and_counts(self) -> None:
        self.assertEqual(self.payload["packet_id"], "P2434")
        self.assertEqual(self.payload["stage_id"], "S1384")
        self.assertEqual(
            self.payload["status"],
            "PASS_CONDITIONAL_LEGACY_ROLE_TRANSFER_CLAIM_LATTICE_NO_ROLE_NO_LTOTAL_NO_TOE_CLOSURE",
        )
        self.assertTrue(OUT.exists())
        self.assertTrue(MD.exists())
        self.assertEqual(self.theorem["role_obligation_count"], 6)
        self.assertEqual(self.theorem["legacy_role_claim_count"], 4)
        self.assertEqual(self.theorem["role_lattice_assignment_count"], 64)

    def test_current_and_minimal_role_lattice(self) -> None:
        self.assertEqual(self.theorem["current_mask"], 0)
        self.assertEqual(self.theorem["current_ready_role_claims"], [])
        self.assertEqual(self.theorem["current_ready_role_claim_count"], 0)
        self.assertEqual(self.theorem["all_roles_ready_masks"], [63])
        self.assertEqual(self.theorem["all_roles_ready_mask_count"], 1)
        self.assertEqual(
            self.theorem["ready_role_count_distribution"],
            {"0": 53, "1": 6, "2": 2, "3": 2, "4": 1},
        )
        self.assertEqual(self.theorem["minimal_size_by_role_claim"]["legacy_weak_mixing_angle"], 3)
        for role_id in [
            "legacy_inverse_alpha_em",
            "legacy_beta_power_gravity_hierarchy",
            "legacy_torsion_to_chi11_orientation",
        ]:
            self.assertEqual(self.theorem["minimal_size_by_role_claim"][role_id], 4)

    def test_necessity_not_sufficiency_and_inheritance(self) -> None:
        self.assertTrue(self.theorem["role_transfer_and_ltotal_are_necessary_for_every_role"])
        self.assertTrue(self.theorem["role_transfer_and_ltotal_are_not_sufficient_for_any_role"])
        self.assertTrue(self.theorem["p2433_convergence_ready_inherited"])
        self.assertTrue(self.theorem["p2433_role_transfer_next_inherited"])

    def test_hard_limits_and_docs(self) -> None:
        self.assertFalse(self.theorem["source_obligation_discharge_exported_by_this_certificate"])
        self.assertFalse(self.theorem["chi11_source_exported_by_this_certificate"])
        self.assertFalse(self.theorem["qw2191_discharged_by_this_certificate"])
        self.assertFalse(self.theorem["role_transfer_licensed"])
        self.assertFalse(self.theorem["role_bearing_ltotal_exported"])
        self.assertFalse(self.theorem["legacy_role_claim_transferred_by_this_certificate"])
        self.assertFalse(self.theorem["toe_closure_exported"])
        self.assertTrue(all(self.payload["gatekeeper_checks"].values()))
        self.assertIn("conditional legacy role-transfer", MD.read_text(encoding="utf-8"))
        self.assertIn("P2434/S1384", (ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md").read_text(encoding="utf-8"))
        self.assertIn("P2434/S1384", (ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md").read_text(encoding="utf-8"))


if __name__ == "__main__":
    unittest.main()
