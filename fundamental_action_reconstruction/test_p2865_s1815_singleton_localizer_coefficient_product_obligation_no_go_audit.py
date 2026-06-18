import json
import subprocess
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
P2864_SCRIPT = ROOT / "p2864_s1814_aut_z12_invariant_nonlocal_boundary_localization_no_go_audit.py"
SCRIPT = ROOT / "p2865_s1815_singleton_localizer_coefficient_product_obligation_no_go_audit.py"
JSON_PATH = ROOT / "generated" / "p2865_s1815_singleton_localizer_coefficient_product_obligation_no_go_audit.json"
MD_PATH = ROOT / "generated" / "p2865_s1815_singleton_localizer_coefficient_product_obligation_no_go_audit.md"


class P2865SingletonLocalizerCoefficientProductObligationNoGoAuditTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        for script in (P2864_SCRIPT, SCRIPT):
            subprocess.run([sys.executable, str(script)], check=True, cwd=ROOT, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
        cls.payload = json.loads(JSON_PATH.read_text(encoding="utf-8"))
        cls.audit = cls.payload["singleton_localizer_coefficient_product_obligation_no_go_audit"]

    def test_status_and_p2864_input(self):
        self.assertEqual(
            self.payload["status"],
            "P2865_SINGLETON_LOCALIZER_COEFFICIENT_PRODUCT_OBLIGATION_NO_GO_AUDIT_NO_CLOSURE",
        )
        self.assertEqual(
            self.audit["input_status_rechecked"],
            "P2864_AUT_Z12_INVARIANT_NONLOCAL_BOUNDARY_LOCALIZATION_NO_GO_AUDIT_NO_CLOSURE",
        )

    def test_unique_exact_singleton_is_d11_with_9_over_5(self):
        exact = self.audit["exact_singleton_rows"]
        self.assertEqual(len(exact), 1)
        self.assertEqual(exact[0]["singleton_d"], 11)
        self.assertEqual(exact[0]["required_scalar"]["fraction"], "9/5")
        self.assertTrue(exact[0]["requires_prime5_coefficient"])
        self.assertTrue(exact[0]["requires_nonpremise_singleton_localizer"])

    def test_all_other_singletons_fail(self):
        failures = [row for row in self.audit["singleton_rows"] if row["singleton_d"] != 11]
        self.assertEqual(len(failures), 10)
        self.assertTrue(all(not row["exact_target_possible_with_scalar"] for row in failures))

    def test_product_obligation_matrix_exports_no_source(self):
        rows = {row["candidate"]: row for row in self.audit["product_obligation_matrix"]}
        self.assertFalse(rows["localizer_without_coefficient_law"]["exports_boundary_source_law"])
        self.assertFalse(rows["coefficient_without_localizer_law"]["exports_boundary_source_law"])
        self.assertTrue(rows["singleton_localizer_times_prime5_coefficient"]["finite_witness_passes"])
        self.assertFalse(rows["singleton_localizer_times_prime5_coefficient"]["exports_boundary_source_law"])
        self.assertFalse(self.payload["acceptance_matrix"]["exports_selector_or_localizer_source"])
        self.assertFalse(self.payload["acceptance_matrix"]["exports_prime5_source_law"])

    def test_no_false_closure_documents(self):
        flags = self.payload["decision"]["negative_export_flags"]
        self.assertFalse(flags["selector_or_localizer_source_exported"])
        self.assertFalse(flags["prime5_source_exported"])
        self.assertFalse(flags["boundary_source_law_exported"])
        self.assertFalse(flags["toe_closure_exported"])
        self.assertIn("P2865/S1815", MD_PATH.read_text(encoding="utf-8"))
        self.assertIn("P2865/S1815", (ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md").read_text(encoding="utf-8"))
        self.assertIn("P2865/S1815", (ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md").read_text(encoding="utf-8"))
        self.assertIn("P2865", (REPO / "AGENTS.md").read_text(encoding="utf-8"))


if __name__ == "__main__":
    unittest.main()
