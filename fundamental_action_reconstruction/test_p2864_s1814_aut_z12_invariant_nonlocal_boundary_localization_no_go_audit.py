import json
import subprocess
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
P2863_SCRIPT = ROOT / "p2863_s1813_dirichlet_boundary_datum_source_class_no_go_audit.py"
SCRIPT = ROOT / "p2864_s1814_aut_z12_invariant_nonlocal_boundary_localization_no_go_audit.py"
JSON_PATH = ROOT / "generated" / "p2864_s1814_aut_z12_invariant_nonlocal_boundary_localization_no_go_audit.json"
MD_PATH = ROOT / "generated" / "p2864_s1814_aut_z12_invariant_nonlocal_boundary_localization_no_go_audit.md"


class P2864AutZ12InvariantNonlocalBoundaryLocalizationNoGoAuditTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        for script in (P2863_SCRIPT, SCRIPT):
            subprocess.run([sys.executable, str(script)], check=True, cwd=ROOT, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
        cls.payload = json.loads(JSON_PATH.read_text(encoding="utf-8"))
        cls.audit = cls.payload["aut_z12_invariant_nonlocal_boundary_localization_no_go_audit"]
        cls.system = cls.audit["linear_system_certificate"]

    def test_status_and_p2863_input(self):
        self.assertEqual(
            self.payload["status"],
            "P2864_AUT_Z12_INVARIANT_NONLOCAL_BOUNDARY_LOCALIZATION_NO_GO_AUDIT_NO_CLOSURE",
        )
        self.assertEqual(
            self.audit["input_status_rechecked"],
            "P2863_DIRICHLET_BOUNDARY_DATUM_SOURCE_CLASS_NO_GO_AUDIT_NO_CLOSURE",
        )

    def test_aut_orbits_and_linear_system_inconsistency(self):
        self.assertEqual(self.system["aut_orbits"], [[1, 5, 7, 11], [2, 10], [3, 9], [4, 8], [6]])
        self.assertEqual(self.system["unit_orbit"], [1, 5, 7, 11])
        self.assertFalse(self.system["linear_system_consistent"])
        self.assertLess(self.system["rank_coefficient_matrix"], self.system["rank_augmented_matrix"])

    def test_unit_orbit_ties_unwanted_prime_components(self):
        self.assertEqual(self.system["unit_orbit_prime5_coefficient"]["fraction"], "1/1")
        self.assertEqual(self.system["unit_orbit_prime7_coefficient"]["fraction"], "1/1")
        self.assertEqual(self.system["unit_orbit_prime11_coefficient"]["fraction"], "1/1")
        self.assertEqual(self.audit["target_prime_exponent_vector"]["11"]["fraction"], "9/5")
        self.assertEqual(self.audit["target_prime_exponent_vector"]["5"]["fraction"], "0/1")
        self.assertEqual(self.audit["target_prime_exponent_vector"]["7"]["fraction"], "0/1")

    def test_singleton_localizer_requires_imports(self):
        singleton = self.audit["selector_localized_singleton_certificate"]
        self.assertTrue(singleton["represents_target_if_both_imported"])
        self.assertTrue(singleton["requires_nonpremise_selector_or_localizer"])
        self.assertTrue(singleton["requires_prime5_coefficient"])
        self.assertFalse(singleton["exports_boundary_source_law"])

    def test_candidate_matrix_and_documents_export_no_closure(self):
        rows = {row["candidate"]: row for row in self.audit["candidate_matrix"]}
        self.assertTrue(rows["aut_z12_invariant_nonlocal_log_localization"]["finite_witness_passes"])
        self.assertFalse(rows["aut_z12_invariant_nonlocal_log_localization"]["exports_boundary_source_law"])
        self.assertTrue(rows["selector_localized_singleton_d11"]["finite_witness_passes"])
        self.assertFalse(rows["selector_localized_singleton_d11"]["exports_eta_source_law"])
        self.assertFalse(self.payload["acceptance_matrix"]["exports_selector_closure"])
        flags = self.payload["decision"]["negative_export_flags"]
        self.assertFalse(flags["selector_closure_exported"])
        self.assertFalse(flags["toe_closure_exported"])
        self.assertIn("P2864/S1814", MD_PATH.read_text(encoding="utf-8"))
        self.assertIn("P2864/S1814", (ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md").read_text(encoding="utf-8"))
        self.assertIn("P2864/S1814", (ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md").read_text(encoding="utf-8"))
        self.assertIn("P2864", (REPO / "AGENTS.md").read_text(encoding="utf-8"))


if __name__ == "__main__":
    unittest.main()
