import json
import subprocess
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
SCRIPT = ROOT / "p2926_s1876_prime_log_value_source_solution_space_certificate.py"
OUT = ROOT / "generated" / "p2926_s1876_prime_log_value_source_solution_space_certificate.json"
MD = ROOT / "generated" / "p2926_s1876_prime_log_value_source_solution_space_certificate.md"


class P2926PrimeLogValueSourceSolutionSpaceCertificateTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        subprocess.run([sys.executable, str(SCRIPT)], cwd=ROOT, check=True, stdout=subprocess.PIPE, text=True)
        cls.payload = json.loads(OUT.read_text(encoding="utf-8"))

    def test_status_and_inputs(self):
        self.assertEqual(self.payload["status"], "P2926_PRIME_LOG_VALUE_SOURCE_SOLUTION_SPACE_CERTIFICATE_NO_ACCEPTED_VALUES")
        self.assertIsNotNone(self.payload["input_hashes"]["P2923"])
        self.assertIsNotNone(self.payload["input_hashes"]["P2925"])
        self.assertTrue(self.payload["acceptance_matrix"]["p2923_factorization_readiness_inherited"])
        self.assertTrue(self.payload["acceptance_matrix"]["p2925_no_delta_unlock_inherited"])

    def test_exact_solution_space_certificate(self):
        cert = self.payload["linear_algebra_certificate"]
        self.assertEqual(cert["variable_count_y1_to_y11"], 11)
        self.assertEqual(cert["additive_equation_count"], 29)
        self.assertEqual(cert["rank_of_additive_character_system"], 6)
        self.assertEqual(cert["nullity_of_additive_character_system"], 5)
        self.assertEqual(cert["free_prime_coordinate_count"], 5)
        self.assertTrue(cert["solution_space_is_prime_value_space"])
        self.assertFalse(cert["prime_values_sourced_by_additivity"])

    def test_candidates_and_no_closure(self):
        acc = self.payload["acceptance_matrix"]
        self.assertEqual(acc["candidate_prime_value_source_count"], 6)
        self.assertEqual(acc["accepted_prime_value_source_count"], 0)
        self.assertTrue(acc["no_new_live_frontier_certificate_exported"])
        self.assertFalse(acc["strict_prime_log_value_source_exported"])
        self.assertFalse(acc["strict_delta_source_law_exported"])
        self.assertFalse(acc["strict_damping_beta_eta_source_exported"])

    def test_false_closure_exports_and_docs(self):
        flags = self.payload["decision"]["negative_export_flags"]
        for key in ["strict_prime_log_value_source_exported", "strict_delta_source_law_exported", "strict_damping_beta_eta_source_exported", "nonproxy_ltotal_exported", "eom_hamiltonian_exported", "bridge_closure_exported", "role_transfer_exported", "toe_closure_exported"]:
            self.assertFalse(flags[key])
        self.assertIn("P2926/S1876", MD.read_text(encoding="utf-8"))
        self.assertIn("P2926/S1876", (ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md").read_text(encoding="utf-8"))
        self.assertIn("P2926/S1876", (ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md").read_text(encoding="utf-8"))
        self.assertIn("P2926", (REPO / "AGENTS.md").read_text(encoding="utf-8"))


if __name__ == "__main__":
    unittest.main()
