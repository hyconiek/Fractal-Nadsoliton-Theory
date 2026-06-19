import json
import subprocess
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
SCRIPT = ROOT / "p2921_s1871_sigma_gamma_unlock_verifier.py"
OUT = ROOT / "generated" / "p2921_s1871_sigma_gamma_unlock_verifier.json"
MD = ROOT / "generated" / "p2921_s1871_sigma_gamma_unlock_verifier.md"


class P2921SigmaGammaUnlockVerifierTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        subprocess.run([sys.executable, str(SCRIPT)], cwd=ROOT, check=True, stdout=subprocess.PIPE, text=True)
        cls.payload = json.loads(OUT.read_text(encoding="utf-8"))

    def test_status_and_p2920_input(self):
        self.assertEqual(self.payload["status"], "P2921_SIGMA_GAMMA_UNLOCK_VERIFIER_EXECUTED_NO_ACCEPTED_SOURCE")
        self.assertTrue(self.payload["acceptance_matrix"]["p2920_rechecked_certificate"])
        self.assertIsNotNone(self.payload["input_hashes"]["P2920"])

    def test_verifier_schema(self):
        schema = self.payload["constructed_theoretical_objects"]["verifier_schema"]
        self.assertEqual(schema["verifier_name"], "Strict_sigma_Gamma_Action_Source_Map_Unlock_Verifier")
        self.assertEqual(len(schema["acceptance_obligations"]), 5)
        self.assertIn("candidate sigma_Gamma", schema["input_type"])

    def test_candidate_counts_and_obligation_pass_counts(self):
        acc = self.payload["acceptance_matrix"]
        self.assertEqual(acc["obligation_count"], 5)
        self.assertEqual(acc["candidate_count"], 6)
        self.assertEqual(acc["accepted_candidate_count"], 0)
        self.assertEqual(acc["computed_nonzero_value_pass_count"], 4)
        self.assertEqual(acc["strict_nadsoliton_provenance_pass_count"], 0)
        self.assertEqual(acc["scale_orbit_breaking_pass_count"], 0)
        self.assertEqual(acc["explicit_IQ_coupling_pass_count"], 6)

    def test_verifier_exported_but_no_closure(self):
        acc = self.payload["acceptance_matrix"]
        self.assertTrue(acc["strict_sigma_gamma_unlock_verifier_exported"])
        self.assertFalse(acc["strict_sigma_gamma_source_accepted"])
        self.assertFalse(acc["accepted_as_nonproxy_ltotal_source"])
        flags = self.payload["decision"]["negative_export_flags"]
        for key in ["nonproxy_ltotal_exported", "eom_closure_exported", "hamiltonian_closure_exported", "bridge_closure_exported", "role_transfer_exported", "toe_closure_exported"]:
            self.assertFalse(flags[key])

    def test_documents_updated(self):
        self.assertIn("P2921/S1871", MD.read_text(encoding="utf-8"))
        self.assertIn("P2921/S1871", (ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md").read_text(encoding="utf-8"))
        self.assertIn("P2921/S1871", (ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md").read_text(encoding="utf-8"))
        self.assertIn("P2921", (REPO / "AGENTS.md").read_text(encoding="utf-8"))


if __name__ == "__main__":
    unittest.main()
