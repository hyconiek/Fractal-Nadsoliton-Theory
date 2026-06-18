import json
import subprocess
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
P2861_SCRIPT = ROOT / "p2861_s1811_log_scale_affine_variational_eta_blind_no_source_audit.py"
SCRIPT = ROOT / "p2862_s1812_log_scale_dirichlet_boundary_eta_recovery_no_source_audit.py"
JSON_PATH = ROOT / "generated" / "p2862_s1812_log_scale_dirichlet_boundary_eta_recovery_no_source_audit.json"
MD_PATH = ROOT / "generated" / "p2862_s1812_log_scale_dirichlet_boundary_eta_recovery_no_source_audit.md"


class P2862LogScaleDirichletBoundaryEtaRecoveryNoSourceAuditTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        for script in (P2861_SCRIPT, SCRIPT):
            subprocess.run([sys.executable, str(script)], check=True, cwd=ROOT, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
        cls.payload = json.loads(JSON_PATH.read_text(encoding="utf-8"))
        cls.audit = cls.payload["log_scale_dirichlet_boundary_eta_recovery_no_source_audit"]

    def test_status_and_p2861_input(self):
        self.assertEqual(
            self.payload["status"],
            "P2862_LOG_SCALE_DIRICHLET_BOUNDARY_ETA_RECOVERY_NO_SOURCE_AUDIT_NO_CLOSURE",
        )
        self.assertEqual(
            self.audit["input_status_rechecked"],
            "P2861_LOG_SCALE_AFFINE_VARIATIONAL_ETA_BLIND_NO_SOURCE_AUDIT_NO_CLOSURE",
        )

    def test_strict_dirichlet_data_recovers_eta_and_beta(self):
        self.assertAlmostEqual(self.audit["recovered_eta_from_strict_boundaries"], 9 / 5)
        self.assertAlmostEqual(self.audit["recovered_beta_from_left_boundary"], 1.0)
        self.assertLess(self.audit["max_reconstruction_residual"], 1e-12)
        self.assertTrue(self.payload["acceptance_matrix"]["facts"]["strict_dirichlet_data_recovers_eta"])
        self.assertTrue(self.payload["acceptance_matrix"]["facts"]["strict_dirichlet_data_recovers_beta"])

    def test_alternative_endpoint_data_selects_other_eta_values(self):
        self.assertGreater(self.audit["non_strict_boundary_eta_count"], 0)
        non_strict_rows = [
            row for row in self.audit["alternative_endpoint_rows"] if not row["is_strict_eta"] and row["abs_eta_error"] < 1e-12
        ]
        self.assertGreater(len(non_strict_rows), 0)
        recovered = {round(row["recovered_eta"], 12) for row in non_strict_rows}
        self.assertIn(1.0, recovered)
        self.assertIn(2.0, recovered)

    def test_candidate_matrix_exports_no_source(self):
        rows = {row["candidate"]: row for row in self.audit["candidate_matrix"]}
        self.assertTrue(rows["strict_endpoint_dirichlet_eta_recovery"]["finite_witness_passes"])
        self.assertFalse(rows["strict_endpoint_dirichlet_eta_recovery"]["exports_eta_source_law"])
        self.assertTrue(rows["alternative_endpoint_counterfamily"]["finite_witness_passes"])
        self.assertFalse(rows["alternative_endpoint_counterfamily"]["exports_boundary_source_law"])
        self.assertEqual(self.audit["accepted_candidate_count"], 0)
        self.assertFalse(self.payload["acceptance_matrix"]["exports_eta_source_law"])
        self.assertFalse(self.payload["acceptance_matrix"]["exports_unit_bearing_coupling_localization_theorem"])

    def test_no_false_closure_documents(self):
        flags = self.payload["decision"]["negative_export_flags"]
        self.assertFalse(flags["eta_source_exported"])
        self.assertFalse(flags["boundary_source_law_exported"])
        self.assertFalse(flags["unit_bearing_coupling_localization_theorem_exported"])
        self.assertFalse(flags["strict_damping_compression_bridge_exported"])
        self.assertFalse(flags["toe_closure_exported"])
        self.assertIn("P2862/S1812", MD_PATH.read_text(encoding="utf-8"))
        self.assertIn("P2862/S1812", (ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md").read_text(encoding="utf-8"))
        self.assertIn("P2862/S1812", (ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md").read_text(encoding="utf-8"))
        self.assertIn("P2862", (REPO / "AGENTS.md").read_text(encoding="utf-8"))


if __name__ == "__main__":
    unittest.main()
