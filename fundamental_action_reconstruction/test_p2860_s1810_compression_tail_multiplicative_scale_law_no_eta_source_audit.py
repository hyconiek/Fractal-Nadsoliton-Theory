import json
import subprocess
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
P2859_SCRIPT = ROOT / "p2859_s1809_eta_beta_profile_identifiability_no_source_audit.py"
SCRIPT = ROOT / "p2860_s1810_compression_tail_multiplicative_scale_law_no_eta_source_audit.py"
JSON_PATH = ROOT / "generated" / "p2860_s1810_compression_tail_multiplicative_scale_law_no_eta_source_audit.json"
MD_PATH = ROOT / "generated" / "p2860_s1810_compression_tail_multiplicative_scale_law_no_eta_source_audit.md"


class P2860CompressionTailMultiplicativeScaleLawNoEtaSourceAuditTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        for script in (P2859_SCRIPT, SCRIPT):
            subprocess.run([sys.executable, str(script)], check=True, cwd=ROOT, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
        cls.payload = json.loads(JSON_PATH.read_text(encoding="utf-8"))
        cls.audit = cls.payload["compression_tail_multiplicative_scale_law_no_eta_source_audit"]

    def test_status_and_p2859_input(self):
        self.assertEqual(
            self.payload["status"],
            "P2860_COMPRESSION_TAIL_MULTIPLICATIVE_SCALE_LAW_NO_ETA_SOURCE_AUDIT_NO_CLOSURE",
        )
        self.assertEqual(
            self.audit["input_status_rechecked"],
            "P2859_ETA_BETA_PROFILE_IDENTIFIABILITY_NO_SOURCE_AUDIT_NO_CLOSURE",
        )

    def test_strict_tail_passes_scale_law(self):
        self.assertGreater(self.audit["product_pair_count"], 0)
        self.assertLess(self.audit["strict_scale_law_max_abs_residual"], 1e-12)
        for row in self.audit["strict_scale_law_rows_first16"]:
            self.assertLess(row["abs_residual"], 1e-12)
        self.assertTrue(self.payload["acceptance_matrix"]["facts"]["strict_tail_passes_multiplicative_scale_law"])

    def test_eta_blindness_and_beta_constraint(self):
        self.assertGreater(self.audit["non_strict_eta_pass_count"], 0)
        non_strict_passes = [
            row for row in self.audit["eta_family_samples"] if row["passes_scale_law"] and not row["is_strict_eta"]
        ]
        self.assertGreater(len(non_strict_passes), 0)
        beta_passes = [
            row for row in self.audit["beta_constraint_samples"] if row["nontrivial_tail"] and row["passes_scale_law"]
        ]
        self.assertEqual([row["beta"]["fraction"] for row in beta_passes], ["1/1"])

    def test_candidate_matrix_exports_no_source(self):
        rows = {row["candidate"]: row for row in self.audit["candidate_matrix"]}
        self.assertTrue(rows["multiplicative_tail_scale_law"]["finite_witness_passes"])
        self.assertFalse(rows["multiplicative_tail_scale_law"]["exports_eta_source_law"])
        self.assertTrue(rows["sampled_eta_continuum_counterexample"]["finite_witness_passes"])
        self.assertFalse(rows["sampled_eta_continuum_counterexample"]["exports_eta_source_law"])
        self.assertEqual(self.audit["accepted_candidate_count"], 0)
        self.assertFalse(self.payload["acceptance_matrix"]["exports_eta_source_law"])
        self.assertFalse(self.payload["acceptance_matrix"]["exports_strict_damping_compression_bridge"])

    def test_no_false_closure_documents(self):
        flags = self.payload["decision"]["negative_export_flags"]
        self.assertFalse(flags["eta_source_exported"])
        self.assertFalse(flags["target_independent_beta_unit_source_exported"])
        self.assertFalse(flags["strict_damping_compression_bridge_exported"])
        self.assertFalse(flags["toe_closure_exported"])
        self.assertIn("P2860/S1810", MD_PATH.read_text(encoding="utf-8"))
        self.assertIn("P2860/S1810", (ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md").read_text(encoding="utf-8"))
        self.assertIn("P2860/S1810", (ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md").read_text(encoding="utf-8"))
        self.assertIn("P2860", (REPO / "AGENTS.md").read_text(encoding="utf-8"))


if __name__ == "__main__":
    unittest.main()
