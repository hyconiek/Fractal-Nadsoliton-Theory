import json
import subprocess
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
SCRIPT = ROOT / "p2927_s1877_strict_damping_beta_eta_source_packet_verifier.py"
OUT = ROOT / "generated" / "p2927_s1877_strict_damping_beta_eta_source_packet_verifier.json"
MD = ROOT / "generated" / "p2927_s1877_strict_damping_beta_eta_source_packet_verifier.md"


class P2927StrictDampingBetaEtaSourcePacketVerifierTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        subprocess.run([sys.executable, str(SCRIPT)], cwd=ROOT, check=True, stdout=subprocess.PIPE, text=True)
        cls.payload = json.loads(OUT.read_text(encoding="utf-8"))

    def test_status_and_inputs(self):
        self.assertEqual(self.payload["status"], "P2927_STRICT_DAMPING_BETA_ETA_SOURCE_PACKET_VERIFIER_NO_ACCEPTED_PACKET")
        self.assertIsNotNone(self.payload["input_hashes"]["P2925"])
        self.assertIsNotNone(self.payload["input_hashes"]["P2926"])
        self.assertTrue(self.payload["acceptance_matrix"]["p2925_delta_source_absent_inherited"])
        self.assertTrue(self.payload["acceptance_matrix"]["p2926_prime_value_source_absent_inherited"])

    def test_finite_verifier_table(self):
        cert = self.payload["finite_verifier_certificate"]
        self.assertEqual(cert["obligation_count"], 4)
        self.assertEqual(cert["status_table_row_count"], 16)
        self.assertEqual(cert["accepting_row_count"], 1)
        self.assertFalse(cert["current_artifact_packet_accepted"])
        self.assertEqual(cert["candidate_packet_count"], 5)
        self.assertEqual(cert["accepted_candidate_packet_count"], 0)

    def test_current_row_and_exports(self):
        current = self.payload["constructed_theoretical_objects"]["current_artifact_row"]
        self.assertFalse(current["strict_prime_log_value_source_Lp"])
        self.assertFalse(current["strict_delta_4_5_eta_9_5_source_law"])
        self.assertFalse(current["strict_beta_eta_coupling_theorem"])
        self.assertTrue(current["strict_nonpromotion_boundary_audit"])
        self.assertFalse(current["accepted_as_strict_damping_beta_eta_source_packet"])
        acc = self.payload["acceptance_matrix"]
        self.assertFalse(acc["strict_damping_beta_eta_source_packet_exported"])

    def test_false_closure_exports_and_docs(self):
        flags = self.payload["decision"]["negative_export_flags"]
        for key in ["strict_prime_log_value_source_exported", "strict_delta_source_law_exported", "strict_beta_eta_coupling_theorem_exported", "strict_damping_beta_eta_source_packet_exported", "nonproxy_ltotal_exported", "eom_hamiltonian_exported", "bridge_closure_exported", "role_transfer_exported", "toe_closure_exported"]:
            self.assertFalse(flags[key])
        self.assertIn("P2927/S1877", MD.read_text(encoding="utf-8"))
        self.assertIn("P2927/S1877", (ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md").read_text(encoding="utf-8"))
        self.assertIn("P2927/S1877", (ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md").read_text(encoding="utf-8"))
        self.assertIn("P2927", (REPO / "AGENTS.md").read_text(encoding="utf-8"))


if __name__ == "__main__":
    unittest.main()
