import json
import subprocess
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
SCRIPT = ROOT / "p2925_s1875_damping_delta_source_linear_system_frontier_certificate.py"
OUT = ROOT / "generated" / "p2925_s1875_damping_delta_source_linear_system_frontier_certificate.json"
MD = ROOT / "generated" / "p2925_s1875_damping_delta_source_linear_system_frontier_certificate.md"


class P2925DampingDeltaSourceLinearSystemFrontierCertificateTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        subprocess.run([sys.executable, str(SCRIPT)], cwd=ROOT, check=True, stdout=subprocess.PIPE, text=True)
        cls.payload = json.loads(OUT.read_text(encoding="utf-8"))

    def test_status_and_inputs(self):
        self.assertEqual(self.payload["status"], "P2925_DAMPING_DELTA_SOURCE_LINEAR_SYSTEM_FRONTIER_CERTIFICATE_NO_NEW_UNLOCK")
        self.assertIsNotNone(self.payload["input_hashes"]["P2923"])
        self.assertIsNotNone(self.payload["input_hashes"]["P2924"])
        self.assertTrue(self.payload["acceptance_matrix"]["p2923_prime_log_readiness_inherited"])
        self.assertTrue(self.payload["acceptance_matrix"]["p2924_no_anchor_inherited"])

    def test_linear_algebra_certificate(self):
        cert = self.payload["linear_algebra_certificate"]
        self.assertEqual(cert["variable_count"], 11)
        self.assertEqual(cert["shape_equation_count"], 10)
        self.assertEqual(cert["rank_without_target_anchor"], 10)
        self.assertEqual(cert["nullity_without_target_anchor"], 1)
        self.assertEqual(cert["rank_with_external_target_anchor"], 11)
        self.assertEqual(cert["nullity_with_external_target_anchor"], 0)
        self.assertTrue(cert["target_anchor_adds_independent_row"])
        self.assertFalse(cert["target_anchor_is_sourced_by_current_rows"])

    def test_minimal_unlock_packet_and_no_closure(self):
        acc = self.payload["acceptance_matrix"]
        self.assertEqual(acc["source_packet_row_count"], 3)
        self.assertEqual(acc["candidate_next_object_count"], 3)
        self.assertEqual(acc["accepted_current_unlock_object_count"], 0)
        self.assertTrue(acc["no_new_live_frontier_certificate_exported"])
        self.assertFalse(acc["strict_prime_log_value_source_exported"])
        self.assertFalse(acc["strict_delta_source_law_exported"])
        self.assertFalse(acc["strict_damping_beta_eta_source_exported"])

    def test_false_closure_exports_and_docs(self):
        flags = self.payload["decision"]["negative_export_flags"]
        for key in ["strict_prime_log_value_source_exported", "strict_delta_source_law_exported", "strict_damping_beta_eta_source_exported", "nonproxy_ltotal_exported", "eom_hamiltonian_exported", "bridge_closure_exported", "role_transfer_exported", "toe_closure_exported"]:
            self.assertFalse(flags[key])
        self.assertIn("P2925/S1875", MD.read_text(encoding="utf-8"))
        self.assertIn("P2925/S1875", (ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md").read_text(encoding="utf-8"))
        self.assertIn("P2925/S1875", (ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md").read_text(encoding="utf-8"))
        self.assertIn("P2925", (REPO / "AGENTS.md").read_text(encoding="utf-8"))


if __name__ == "__main__":
    unittest.main()
