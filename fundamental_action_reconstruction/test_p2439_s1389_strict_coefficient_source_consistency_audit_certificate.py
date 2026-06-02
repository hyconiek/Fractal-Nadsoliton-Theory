import json
import subprocess
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent
SCRIPT = ROOT / "p2439_s1389_strict_coefficient_source_consistency_audit_certificate.py"
OUT = ROOT / "generated" / "p2439_s1389_strict_coefficient_source_consistency_audit_certificate.json"
MD = ROOT / "generated" / "p2439_s1389_strict_coefficient_source_consistency_audit_certificate.md"
P2438 = ROOT / "generated" / "p2438_s1388_strict_kernel_sm_gr_generation_obligation_matrix_certificate.json"


class P2439StrictCoefficientSourceConsistencyAuditCertificateTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        if not P2438.exists():
            subprocess.run([sys.executable, str(ROOT / "p2438_s1388_strict_kernel_sm_gr_generation_obligation_matrix_certificate.py")], check=True)
        subprocess.run([sys.executable, str(SCRIPT)], check=True)
        cls.payload = json.loads(OUT.read_text(encoding="utf-8"))
        cls.cert = cls.payload["strict_coefficient_source_consistency_audit_certificate"]
        cls.theorem = cls.cert["theorem_export"]

    def test_identity_outputs_and_audit_size(self) -> None:
        self.assertEqual(self.payload["packet_id"], "P2439")
        self.assertEqual(self.payload["stage_id"], "S1389")
        self.assertEqual(
            self.payload["status"],
            "PASS_STRICT_COEFFICIENT_SOURCE_AUDIT_NO_CURRENT_FULL_SM_GR_VALUE_GENERATOR",
        )
        self.assertTrue(OUT.exists())
        self.assertTrue(MD.exists())
        self.assertEqual(self.theorem["audited_source_count"], 3)
        self.assertEqual(self.theorem["feature_count"], 9)
        self.assertEqual(len(self.theorem["source_feature_rows"]), 3)
        self.assertEqual(len(self.theorem["source_feature_incidence_matrix_gf2"]), 3)

    def test_rank_tuple_and_no_generator_results(self) -> None:
        self.assertEqual(self.theorem["source_feature_incidence_rank_gf2"], 3)
        self.assertEqual(self.theorem["current_tuple_matching_source_ids"], ["P1563_effective_three_coefficient_chain"])
        self.assertEqual(self.theorem["current_tuple_matching_source_count"], 1)
        self.assertEqual(self.theorem["current_tuple_full_sm_gr_coefficient_source_ids"], [])
        self.assertEqual(self.theorem["current_tuple_full_sm_gr_coefficient_source_count"], 0)
        self.assertEqual(self.theorem["acceptable_current_strict_physical_value_generator_source_ids"], [])
        self.assertEqual(self.theorem["acceptable_current_strict_physical_value_generator_source_count"], 0)

    def test_source_specific_classifications_and_hard_limits(self) -> None:
        self.assertTrue(self.theorem["p1563_current_tuple_but_three_effective_coefficients_only"])
        self.assertEqual(self.theorem["p1664_full_manifest_tuple_mismatch_count"], 3)
        self.assertTrue(self.theorem["p1664_inverse_recovery_is_only_local"])
        self.assertTrue(self.theorem["p1910_all_coefficients_open_symbolic_export"])
        self.assertTrue(self.theorem["p2438_no_strict_generation_inherited"])
        self.assertFalse(self.theorem["strict_kernel_to_coefficient_map_theorem_exported_by_this_certificate"])
        self.assertFalse(self.theorem["strict_observable_value_generator_exported"])
        self.assertFalse(self.theorem["qw2191_discharged"])
        self.assertFalse(self.theorem["role_bearing_ltotal_exported"])
        self.assertFalse(self.theorem["toe_closure_exported"])
        self.assertTrue(all(self.payload["gatekeeper_checks"].values()))

    def test_docs(self) -> None:
        self.assertIn("strict coefficient-source consistency audit", MD.read_text(encoding="utf-8"))
        self.assertIn("P2439/S1389", (ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md").read_text(encoding="utf-8"))
        self.assertIn("P2439/S1389", (ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md").read_text(encoding="utf-8"))


if __name__ == "__main__":
    unittest.main()
