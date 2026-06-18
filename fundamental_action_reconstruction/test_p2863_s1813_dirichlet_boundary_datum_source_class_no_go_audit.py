import json
import subprocess
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
P2862_SCRIPT = ROOT / "p2862_s1812_log_scale_dirichlet_boundary_eta_recovery_no_source_audit.py"
SCRIPT = ROOT / "p2863_s1813_dirichlet_boundary_datum_source_class_no_go_audit.py"
JSON_PATH = ROOT / "generated" / "p2863_s1813_dirichlet_boundary_datum_source_class_no_go_audit.json"
MD_PATH = ROOT / "generated" / "p2863_s1813_dirichlet_boundary_datum_source_class_no_go_audit.md"


class P2863DirichletBoundaryDatumSourceClassNoGoAuditTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        for script in (P2862_SCRIPT, SCRIPT):
            subprocess.run([sys.executable, str(script)], check=True, cwd=ROOT, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
        cls.payload = json.loads(JSON_PATH.read_text(encoding="utf-8"))
        cls.audit = cls.payload["dirichlet_boundary_datum_source_class_no_go_audit"]

    def test_status_and_p2862_input(self):
        self.assertEqual(
            self.payload["status"],
            "P2863_DIRICHLET_BOUNDARY_DATUM_SOURCE_CLASS_NO_GO_AUDIT_NO_CLOSURE",
        )
        self.assertEqual(
            self.audit["input_status_rechecked"],
            "P2862_LOG_SCALE_DIRICHLET_BOUNDARY_ETA_RECOVERY_NO_SOURCE_AUDIT_NO_CLOSURE",
        )

    def test_target_prime_exponent_and_z12_denominator_obstruction(self):
        self.assertEqual(self.audit["target_prime_exponent_vector"]["11"]["fraction"], "9/5")
        obstruction = self.audit["denominator_obstruction"]
        self.assertEqual(obstruction["target_exponent_denominator_support"], [5])
        self.assertEqual(obstruction["allowed_denominator_support"], [2, 3])
        self.assertFalse(obstruction["can_represent_target_exponent_without_prime5"])

    def test_coefficient_and_moment_scans_do_not_source_boundary(self):
        pure_z12, prime5 = self.audit["coefficient_scan_rows"]
        self.assertEqual(pure_z12["exact_match_count"], 0)
        self.assertEqual(prime5["exact_match_count"], 1)
        self.assertTrue(prime5["imports_prime5"])
        self.assertEqual(self.audit["integer_moment_scan"]["exact_match_count"], 0)
        self.assertEqual(self.audit["integer_moment_scan"]["best_diff"]["11"]["fraction"], "1/5")

    def test_candidate_matrix_exports_no_source(self):
        rows = {row["candidate"]: row for row in self.audit["candidate_matrix"]}
        self.assertTrue(rows["pure_z12_smooth_boundary_coefficient"]["finite_witness_passes"])
        self.assertFalse(rows["pure_z12_smooth_boundary_coefficient"]["exports_boundary_source_law"])
        self.assertTrue(rows["prime5_extended_boundary_coefficient"]["finite_witness_passes"])
        self.assertFalse(rows["prime5_extended_boundary_coefficient"]["exports_boundary_source_law"])
        self.assertEqual(self.audit["accepted_candidate_count"], 0)
        self.assertFalse(self.payload["acceptance_matrix"]["exports_boundary_source_law"])
        self.assertFalse(self.payload["acceptance_matrix"]["exports_eta_source_law"])

    def test_no_false_closure_documents(self):
        flags = self.payload["decision"]["negative_export_flags"]
        self.assertFalse(flags["boundary_source_law_exported"])
        self.assertFalse(flags["eta_source_exported"])
        self.assertFalse(flags["prime5_source_exported"])
        self.assertFalse(flags["unit_bearing_coupling_localization_theorem_exported"])
        self.assertFalse(flags["toe_closure_exported"])
        self.assertIn("P2863/S1813", MD_PATH.read_text(encoding="utf-8"))
        self.assertIn("P2863/S1813", (ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md").read_text(encoding="utf-8"))
        self.assertIn("P2863/S1813", (ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md").read_text(encoding="utf-8"))
        self.assertIn("P2863", (REPO / "AGENTS.md").read_text(encoding="utf-8"))


if __name__ == "__main__":
    unittest.main()
