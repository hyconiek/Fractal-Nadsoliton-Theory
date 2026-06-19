import json
import subprocess
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
SCRIPT = ROOT / "p2934_s1884_aut_breaking_prime_coordinate_source_law_acceptance_verifier.py"
OUT = ROOT / "generated" / "p2934_s1884_aut_breaking_prime_coordinate_source_law_acceptance_verifier.json"
MD = ROOT / "generated" / "p2934_s1884_aut_breaking_prime_coordinate_source_law_acceptance_verifier.md"


class P2934AutBreakingPrimeCoordinateSourceLawAcceptanceVerifierTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        subprocess.run([sys.executable, str(SCRIPT)], cwd=ROOT, check=True, stdout=subprocess.PIPE, text=True)
        cls.payload = json.loads(OUT.read_text(encoding="utf-8"))

    def test_status_and_input(self):
        self.assertEqual(self.payload["status"], "P2934_AUT_BREAKING_SOURCE_LAW_ACCEPTANCE_VERIFIER_NO_ACCEPTED_SOURCE")
        self.assertIsNotNone(self.payload["input_hashes"]["P2933"])

    def test_acceptance_certificate(self):
        cert = self.payload["acceptance_certificate"]
        self.assertEqual(cert["obligation_count"], 5)
        self.assertEqual(cert["satisfied_obligation_count"], 3)
        self.assertEqual(cert["bounded_total_vector_count"], 243)
        self.assertEqual(cert["bounded_nonzero_vector_count"], 242)
        self.assertEqual(cert["bounded_product_additive_nonzero_count"], 242)
        self.assertEqual(cert["bounded_aut_breaking_nonzero_count"], 242)
        self.assertEqual(cert["accepted_strict_source_law_count"], 0)
        self.assertTrue(cert["all_current_candidates_rejected"])

    def test_verifier_and_no_closure_flags(self):
        verifier = self.payload["constructed_theoretical_objects"]["verifier"]
        self.assertEqual(verifier["name"], "Strict_AutBreaking_PrimeCoordinate_Source_Law_Acceptance_Verifier")
        self.assertEqual(verifier["target_future_object"], "Strict_AutBreaking_PrimeCoordinate_Source_Law")
        flags = self.payload["decision"]["negative_export_flags"]
        for key in ["strict_prime_log_value_source_exported", "strict_delta_eta_source_exported", "strict_beta_eta_coupling_theorem_exported", "strict_damping_beta_eta_source_packet_exported", "nonproxy_ltotal_exported", "eom_hamiltonian_exported", "bridge_closure_exported", "role_transfer_exported", "toe_closure_exported"]:
            self.assertFalse(flags[key])

    def test_docs_updated(self):
        self.assertIn("P2934/S1884", MD.read_text(encoding="utf-8"))
        self.assertIn("P2934/S1884", (ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md").read_text(encoding="utf-8"))
        self.assertIn("P2934/S1884", (ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md").read_text(encoding="utf-8"))
        self.assertIn("P2934", (REPO / "AGENTS.md").read_text(encoding="utf-8"))


if __name__ == "__main__":
    unittest.main()
