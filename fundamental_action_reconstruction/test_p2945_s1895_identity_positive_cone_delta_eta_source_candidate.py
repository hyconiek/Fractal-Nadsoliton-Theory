import json
import subprocess
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
SCRIPT = ROOT / "p2945_s1895_identity_positive_cone_delta_eta_source_candidate.py"
OUT = ROOT / "generated" / "p2945_s1895_identity_positive_cone_delta_eta_source_candidate.json"
MD = ROOT / "generated" / "p2945_s1895_identity_positive_cone_delta_eta_source_candidate.md"


class P2945IdentityPositiveConeDeltaEtaSourceCandidateTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        subprocess.run([sys.executable, str(SCRIPT)], cwd=ROOT, check=True, stdout=subprocess.PIPE, text=True)
        cls.payload = json.loads(OUT.read_text(encoding="utf-8"))

    def test_status_and_input(self):
        self.assertEqual(self.payload["status"], "P2945_IDENTITY_POSITIVE_CONE_DELTA_ETA_SOURCE_CANDIDATE_CONDITIONAL_ONLY")
        self.assertIsNotNone(self.payload["input_hashes"]["P2944"])

    def test_delta_eta_candidate_certificate(self):
        cert = self.payload["delta_eta_candidate_certificate"]
        self.assertTrue(cert["finite_candidate_passes_arithmetic_gate"])
        self.assertEqual(cert["delta_candidate"]["as_string"], "4/5")
        self.assertEqual(cert["eta_candidate"]["as_string"], "9/5")
        self.assertTrue(cert["eta_equals_one_plus_delta"])
        self.assertFalse(cert["strict_ratio_source_theorem_exported"])
        self.assertFalse(cert["strict_beta_eta_coupling_theorem_exported"])
        self.assertFalse(cert["accepted_strict_delta_eta_source_law"])

    def test_formula_rows_acceptance_and_nonpromotion(self):
        obj = self.payload["constructed_theoretical_objects"]
        inventory = obj["source_inventory"]
        self.assertEqual(inventory["prime_count"], 5)
        self.assertEqual(inventory["identity_count"], 1)
        self.assertEqual(inventory["prime_vector_sum"], 9)
        formulas = obj["candidate_formula_rows"]
        self.assertTrue(formulas[0]["matches_target_delta_4_5"])
        self.assertTrue(formulas[1]["matches_target_eta_9_5"])
        self.assertTrue(formulas[2]["matches_eta"])
        acceptance = obj["acceptance_rows"]
        self.assertTrue(all(row["satisfied"] for row in acceptance[:4]))
        self.assertFalse(acceptance[4]["satisfied"])
        self.assertFalse(acceptance[5]["satisfied"])
        flags = self.payload["decision"]["negative_export_flags"]
        for key in ["strict_identity_grounding_theorem_exported", "strict_ratio_source_theorem_exported", "strict_prime_log_value_source_exported", "strict_delta_eta_source_law_exported", "strict_beta_eta_coupling_theorem_exported", "strict_damping_beta_eta_source_packet_exported", "nonproxy_ltotal_exported", "eom_hamiltonian_exported", "bridge_closure_exported", "role_transfer_exported", "toe_closure_exported"]:
            self.assertFalse(flags[key])

    def test_docs_updated(self):
        self.assertIn("P2945/S1895", MD.read_text(encoding="utf-8"))
        self.assertIn("P2945/S1895", (ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md").read_text(encoding="utf-8"))
        self.assertIn("P2945/S1895", (ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md").read_text(encoding="utf-8"))
        self.assertIn("P2945", (REPO / "AGENTS.md").read_text(encoding="utf-8"))


if __name__ == "__main__":
    unittest.main()
