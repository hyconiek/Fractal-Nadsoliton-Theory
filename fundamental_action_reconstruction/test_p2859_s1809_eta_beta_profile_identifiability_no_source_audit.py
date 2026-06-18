import json
import subprocess
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
P2858_SCRIPT = ROOT / "p2858_s1808_phase_bit_cell_continuum_no_source_audit.py"
SCRIPT = ROOT / "p2859_s1809_eta_beta_profile_identifiability_no_source_audit.py"
JSON_PATH = ROOT / "generated" / "p2859_s1809_eta_beta_profile_identifiability_no_source_audit.json"
MD_PATH = ROOT / "generated" / "p2859_s1809_eta_beta_profile_identifiability_no_source_audit.md"


class P2859EtaBetaProfileIdentifiabilityNoSourceAuditTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        for script in (P2858_SCRIPT, SCRIPT):
            subprocess.run([sys.executable, str(script)], check=True, cwd=ROOT, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
        cls.payload = json.loads(JSON_PATH.read_text(encoding="utf-8"))
        cls.audit = cls.payload["eta_beta_profile_identifiability_no_source_audit"]

    def test_status_and_p2858_input(self):
        self.assertEqual(self.payload["status"], "P2859_ETA_BETA_PROFILE_IDENTIFIABILITY_NO_SOURCE_AUDIT_NO_CLOSURE")
        self.assertEqual(
            self.audit["input_status_rechecked"],
            "P2858_PHASE_BIT_CELL_CONTINUUM_NO_SOURCE_AUDIT_NO_CLOSURE",
        )

    def test_two_point_inverse_recovers_supplied_parameters(self):
        inverse = self.audit["two_point_inverse"]
        self.assertEqual(inverse["d_pair"], [2, 3])
        self.assertLess(inverse["eta_abs_error"], 1e-12)
        self.assertLess(inverse["beta_abs_error"], 1e-12)
        facts = self.payload["acceptance_matrix"]["facts"]
        self.assertTrue(facts["two_point_inverse_recovers_eta"])
        self.assertTrue(facts["two_point_inverse_recovers_beta"])

    def test_jacobian_and_eta_prime5_boundary(self):
        jac = self.audit["jacobian_certificate"]
        self.assertTrue(jac["nonzero_determinant"])
        self.assertTrue(jac["local_identifiability"])
        eta = self.audit["strict_eta"]
        self.assertEqual(eta["eta"]["fraction"], "9/5")
        self.assertTrue(eta["requires_prime5"])
        self.assertFalse(eta["pure_z12_denominator_source"])

    def test_candidate_matrix_exports_no_source(self):
        rows = {row["candidate"]: row for row in self.audit["candidate_matrix"]}
        self.assertTrue(rows["two_point_profile_inverse"]["finite_witness_passes"])
        self.assertFalse(rows["two_point_profile_inverse"]["exports_pre_profile_source_law"])
        self.assertTrue(rows["jacobian_local_identifiability"]["finite_witness_passes"])
        self.assertFalse(rows["jacobian_local_identifiability"]["exports_pre_profile_source_law"])
        self.assertEqual(self.audit["accepted_candidate_count"], 0)
        self.assertFalse(self.payload["acceptance_matrix"]["exports_eta_beta_source_law"])
        self.assertFalse(self.payload["acceptance_matrix"]["exports_strict_damping_compression_bridge"])

    def test_no_false_closure_documents(self):
        flags = self.payload["decision"]["negative_export_flags"]
        self.assertFalse(flags["eta_beta_source_law_exported"])
        self.assertFalse(flags["strict_damping_compression_bridge_exported"])
        self.assertFalse(flags["full_kernel_bridge_exported"])
        self.assertFalse(flags["toe_closure_exported"])
        self.assertIn("P2859/S1809", MD_PATH.read_text(encoding="utf-8"))
        self.assertIn("P2859/S1809", (ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md").read_text(encoding="utf-8"))
        self.assertIn("P2859/S1809", (ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md").read_text(encoding="utf-8"))
        self.assertIn("P2859", (REPO / "AGENTS.md").read_text(encoding="utf-8"))


if __name__ == "__main__":
    unittest.main()
