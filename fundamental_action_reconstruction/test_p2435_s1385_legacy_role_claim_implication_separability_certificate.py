import json
import subprocess
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent
SCRIPT = ROOT / "p2435_s1385_legacy_role_claim_implication_separability_certificate.py"
OUT = ROOT / "generated" / "p2435_s1385_legacy_role_claim_implication_separability_certificate.json"
MD = ROOT / "generated" / "p2435_s1385_legacy_role_claim_implication_separability_certificate.md"
P2434 = ROOT / "generated" / "p2434_s1384_conditional_legacy_role_transfer_claim_lattice_certificate.json"


class P2435LegacyRoleClaimImplicationSeparabilityCertificateTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        if not P2434.exists():
            subprocess.run([sys.executable, str(ROOT / "p2434_s1384_conditional_legacy_role_transfer_claim_lattice_certificate.py")], check=True)
        subprocess.run([sys.executable, str(SCRIPT)], check=True)
        cls.payload = json.loads(OUT.read_text(encoding="utf-8"))
        cls.cert = cls.payload["legacy_role_claim_implication_separability_certificate"]
        cls.theorem = cls.cert["theorem_export"]

    def test_identity_outputs_and_rank(self) -> None:
        self.assertEqual(self.payload["packet_id"], "P2435")
        self.assertEqual(self.payload["stage_id"], "S1385")
        self.assertEqual(
            self.payload["status"],
            "PASS_LEGACY_ROLE_CLAIM_IMPLICATION_SEPARABILITY_NO_ROLE_TRANSFER_NO_CLOSURE",
        )
        self.assertTrue(OUT.exists())
        self.assertTrue(MD.exists())
        self.assertEqual(self.theorem["role_claim_count"], 4)
        self.assertEqual(self.theorem["role_obligation_count"], 6)
        self.assertEqual(self.theorem["incidence_rank_gf2"], 4)

    def test_implication_poset(self) -> None:
        self.assertEqual(self.theorem["implication_row_count"], 16)
        self.assertEqual(self.theorem["implication_true_count"], 5)
        self.assertEqual(
            self.theorem["nontrivial_implications"],
            [["legacy_inverse_alpha_em", "legacy_weak_mixing_angle"]],
        )
        self.assertEqual(self.theorem["nonimplication_witness_count"], 11)
        self.assertTrue(self.theorem["alpha_em_strictly_implies_weak_mixing_obligation_readiness"])
        self.assertTrue(self.theorem["weak_mixing_does_not_imply_alpha_em"])
        self.assertTrue(self.theorem["no_two_distinct_claims_equivalent"])

    def test_coverage_and_hard_limits(self) -> None:
        coverage = self.theorem["obligation_coverage_by_role"]
        self.assertEqual(len(coverage["role_transfer_audit_license"]), 4)
        self.assertEqual(len(coverage["role_bearing_ltotal_export"]), 4)
        self.assertEqual(len(coverage["alpha_geo_strict_role_successor_theorem"]), 2)
        self.assertEqual(len(coverage["beta_tors_strict_role_successor_theorem"]), 3)
        self.assertEqual(len(coverage["strict_nonlinear_hierarchy_successor_theorem"]), 1)
        self.assertEqual(len(coverage["chi11_orientation_role_successor_theorem"]), 1)
        self.assertTrue(self.theorem["role_transfer_and_ltotal_cover_all_claims_but_do_not_separate_them"])
        self.assertFalse(self.theorem["legacy_role_claim_transferred_by_this_certificate"])
        self.assertFalse(self.theorem["role_transfer_licensed"])
        self.assertFalse(self.theorem["role_bearing_ltotal_exported"])
        self.assertFalse(self.theorem["toe_closure_exported"])
        self.assertTrue(all(self.payload["gatekeeper_checks"].values()))

    def test_docs(self) -> None:
        self.assertIn("legacy role-claim implication", MD.read_text(encoding="utf-8"))
        self.assertIn("P2435/S1385", (ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md").read_text(encoding="utf-8"))
        self.assertIn("P2435/S1385", (ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md").read_text(encoding="utf-8"))


if __name__ == "__main__":
    unittest.main()
