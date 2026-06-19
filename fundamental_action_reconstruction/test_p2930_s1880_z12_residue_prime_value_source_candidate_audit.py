import json
import subprocess
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
SCRIPT = ROOT / "p2930_s1880_z12_residue_prime_value_source_candidate_audit.py"
OUT = ROOT / "generated" / "p2930_s1880_z12_residue_prime_value_source_candidate_audit.json"
MD = ROOT / "generated" / "p2930_s1880_z12_residue_prime_value_source_candidate_audit.md"


class P2930Z12ResiduePrimeValueSourceCandidateAuditTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        subprocess.run([sys.executable, str(SCRIPT)], cwd=ROOT, check=True, stdout=subprocess.PIPE, text=True)
        cls.payload = json.loads(OUT.read_text(encoding="utf-8"))

    def test_status_and_input(self):
        self.assertEqual(self.payload["status"], "P2930_Z12_RESIDUE_PRIME_VALUE_SOURCE_CANDIDATE_AUDIT_REJECTED_AS_STRICT_SOURCE")
        self.assertIsNotNone(self.payload["input_hashes"]["P2929"])

    def test_candidate_values_and_finite_carrier(self):
        candidate = self.payload["constructed_theoretical_objects"]["candidate"]
        self.assertEqual(candidate["prime_values"], {"L_2": 2, "L_3": 3, "L_5": 5, "L_7": 5, "L_11": 1})
        matrix = self.payload["acceptance_matrix"]
        self.assertTrue(matrix["computes_five_nonzero_prime_labels"])
        self.assertEqual(matrix["audited_product_pair_count"], 29)
        self.assertEqual(matrix["formal_additive_defect_count"], 0)
        self.assertFalse(matrix["accepted_as_strict_prime_log_value_source"])

    def test_rejection_and_no_closure_flags(self):
        matrix = self.payload["acceptance_matrix"]
        self.assertFalse(matrix["strict_nadsoliton_source_theorem_exported"])
        self.assertFalse(matrix["intrinsic_nonconventional_scale_exported"])
        flags = self.payload["decision"]["negative_export_flags"]
        for key in ["strict_prime_log_value_source_exported", "strict_delta_eta_source_exported", "strict_beta_eta_coupling_theorem_exported", "strict_damping_beta_eta_source_packet_exported", "nonproxy_ltotal_exported", "eom_hamiltonian_exported", "bridge_closure_exported", "role_transfer_exported", "toe_closure_exported"]:
            self.assertFalse(flags[key])

    def test_docs_updated(self):
        self.assertIn("P2930/S1880", MD.read_text(encoding="utf-8"))
        self.assertIn("P2930/S1880", (ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md").read_text(encoding="utf-8"))
        self.assertIn("P2930/S1880", (ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md").read_text(encoding="utf-8"))
        self.assertIn("P2930", (REPO / "AGENTS.md").read_text(encoding="utf-8"))


if __name__ == "__main__":
    unittest.main()
