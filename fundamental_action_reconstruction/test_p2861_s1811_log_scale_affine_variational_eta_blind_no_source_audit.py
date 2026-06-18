import json
import subprocess
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
P2860_SCRIPT = ROOT / "p2860_s1810_compression_tail_multiplicative_scale_law_no_eta_source_audit.py"
SCRIPT = ROOT / "p2861_s1811_log_scale_affine_variational_eta_blind_no_source_audit.py"
JSON_PATH = ROOT / "generated" / "p2861_s1811_log_scale_affine_variational_eta_blind_no_source_audit.json"
MD_PATH = ROOT / "generated" / "p2861_s1811_log_scale_affine_variational_eta_blind_no_source_audit.md"


class P2861LogScaleAffineVariationalEtaBlindNoSourceAuditTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        for script in (P2860_SCRIPT, SCRIPT):
            subprocess.run([sys.executable, str(script)], check=True, cwd=ROOT, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
        cls.payload = json.loads(JSON_PATH.read_text(encoding="utf-8"))
        cls.audit = cls.payload["log_scale_affine_variational_eta_blind_no_source_audit"]

    def test_status_and_p2860_input(self):
        self.assertEqual(
            self.payload["status"],
            "P2861_LOG_SCALE_AFFINE_VARIATIONAL_ETA_BLIND_NO_SOURCE_AUDIT_NO_CLOSURE",
        )
        self.assertEqual(
            self.audit["input_status_rechecked"],
            "P2860_COMPRESSION_TAIL_MULTIPLICATIVE_SCALE_LAW_NO_ETA_SOURCE_AUDIT_NO_CLOSURE",
        )

    def test_strict_tail_has_zero_log_scale_affine_action(self):
        self.assertGreater(self.audit["triple_count"], 0)
        self.assertLess(self.audit["strict_scale_affine_action"], 1e-12)
        for row in self.audit["strict_residual_rows_first16"]:
            self.assertLess(row["abs_residual"], 1e-12)
        self.assertTrue(
            self.payload["acceptance_matrix"]["facts"]["strict_tail_satisfies_log_scale_affine_variational_equation"]
        )

    def test_eta_and_beta_blindness(self):
        self.assertGreater(self.audit["non_strict_eta_pass_count"], 0)
        self.assertGreater(self.audit["non_strict_beta_pass_count"], 0)
        non_strict_eta_passes = [
            row for row in self.audit["eta_action_samples"] if row["passes_variational_equation"] and not row["is_strict_eta"]
        ]
        non_strict_beta_passes = [
            row for row in self.audit["beta_action_samples"] if row["passes_variational_equation"] and not row["is_strict_beta"]
        ]
        self.assertGreater(len(non_strict_eta_passes), 0)
        self.assertGreater(len(non_strict_beta_passes), 0)

    def test_candidate_matrix_exports_no_source(self):
        rows = {row["candidate"]: row for row in self.audit["candidate_matrix"]}
        self.assertTrue(rows["log_scale_affine_variational_equation"]["finite_witness_passes"])
        self.assertFalse(rows["log_scale_affine_variational_equation"]["exports_eta_source_law"])
        self.assertTrue(rows["eta_selection_by_scale_affine_action"]["finite_witness_passes"])
        self.assertFalse(rows["eta_selection_by_scale_affine_action"]["exports_eta_source_law"])
        self.assertEqual(self.audit["accepted_candidate_count"], 0)
        self.assertFalse(self.payload["acceptance_matrix"]["exports_eta_source_law"])
        self.assertFalse(self.payload["acceptance_matrix"]["exports_unit_bearing_variational_source_theorem"])

    def test_no_false_closure_documents(self):
        flags = self.payload["decision"]["negative_export_flags"]
        self.assertFalse(flags["eta_source_exported"])
        self.assertFalse(flags["target_independent_beta_unit_source_exported"])
        self.assertFalse(flags["unit_bearing_variational_source_theorem_exported"])
        self.assertFalse(flags["strict_damping_compression_bridge_exported"])
        self.assertFalse(flags["toe_closure_exported"])
        self.assertIn("P2861/S1811", MD_PATH.read_text(encoding="utf-8"))
        self.assertIn("P2861/S1811", (ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md").read_text(encoding="utf-8"))
        self.assertIn("P2861/S1811", (ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md").read_text(encoding="utf-8"))
        self.assertIn("P2861", (REPO / "AGENTS.md").read_text(encoding="utf-8"))


if __name__ == "__main__":
    unittest.main()
