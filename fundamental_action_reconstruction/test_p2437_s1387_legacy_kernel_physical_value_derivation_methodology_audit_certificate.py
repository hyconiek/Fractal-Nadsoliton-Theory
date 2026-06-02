import json
import math
import subprocess
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent
SCRIPT = ROOT / "p2437_s1387_legacy_kernel_physical_value_derivation_methodology_audit_certificate.py"
OUT = ROOT / "generated" / "p2437_s1387_legacy_kernel_physical_value_derivation_methodology_audit_certificate.json"
MD = ROOT / "generated" / "p2437_s1387_legacy_kernel_physical_value_derivation_methodology_audit_certificate.md"
P2436 = ROOT / "generated" / "p2436_s1386_claim_specific_successor_frontier_antichain_certificate.json"


class P2437LegacyKernelPhysicalValueMethodologyAuditCertificateTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        if not P2436.exists():
            subprocess.run([sys.executable, str(ROOT / "p2436_s1386_claim_specific_successor_frontier_antichain_certificate.py")], check=True)
        subprocess.run([sys.executable, str(SCRIPT)], check=True)
        cls.payload = json.loads(OUT.read_text(encoding="utf-8"))
        cls.cert = cls.payload["legacy_kernel_physical_value_derivation_methodology_audit_certificate"]
        cls.theorem = cls.cert["theorem_export"]

    def test_identity_outputs_and_numeric_formulas(self) -> None:
        self.assertEqual(self.payload["packet_id"], "P2437")
        self.assertEqual(self.payload["stage_id"], "S1387")
        self.assertEqual(
            self.payload["status"],
            "PASS_LEGACY_KERNEL_PHYSICAL_VALUE_METHODOLOGY_AUDIT_NO_STRICT_VALUE_GENERATOR_NO_BETA_TORS_CHI11_THEOREM",
        )
        self.assertTrue(OUT.exists())
        self.assertTrue(MD.exists())
        alpha_geo = 4.0 * math.log(2.0)
        beta_tors = 0.01
        self.assertAlmostEqual(self.theorem["legacy_weak_mixing_numeric"], alpha_geo / 12.0)
        self.assertAlmostEqual(self.theorem["legacy_inverse_alpha_em_numeric"], alpha_geo / (2.0 * beta_tors) * (1.0 - beta_tors))
        self.assertAlmostEqual(self.theorem["legacy_beta_power_n20_numeric"], beta_tors**20)

    def test_methodology_reclassification(self) -> None:
        self.assertEqual(self.theorem["audited_legacy_claim_count"], 4)
        self.assertTrue(self.theorem["all_legacy_claims_reclassified_not_strict_generated"])
        self.assertTrue(self.theorem["all_legacy_claims_not_physically_sufficient_from_legacy_kernel_only"])
        self.assertTrue(self.theorem["beta_tors_to_chi11_reclassified_as_search_assumption"])
        self.assertFalse(self.theorem["strict_kernel_physical_value_generator_exported_in_current_repo"])
        self.assertTrue(self.theorem["strict_kernel_should_be_generation_source_for_physical_values"])
        self.assertTrue(self.theorem["k1_kernel_split_guardrail_detected"])
        self.assertTrue(self.theorem["s2_no_silent_role_transfer_detected"])
        self.assertTrue(self.theorem["p2436_claim_specific_frontier_inherited"])

    def test_hard_limits_and_docs(self) -> None:
        self.assertFalse(self.theorem["legacy_role_claim_transferred_by_this_certificate"])
        self.assertFalse(self.theorem["beta_tors_chi11_theorem_exported"])
        self.assertFalse(self.theorem["strict_physical_value_theorem_exported"])
        self.assertFalse(self.theorem["role_bearing_ltotal_exported"])
        self.assertFalse(self.theorem["toe_closure_exported"])
        self.assertTrue(all(self.payload["gatekeeper_checks"].values()))
        self.assertIn("legacy-kernel physical-value", MD.read_text(encoding="utf-8"))
        self.assertIn("P2437/S1387", (ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md").read_text(encoding="utf-8"))
        self.assertIn("P2437/S1387", (ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md").read_text(encoding="utf-8"))


if __name__ == "__main__":
    unittest.main()
