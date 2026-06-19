import json
import subprocess
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
SCRIPT = ROOT / "p2922_s1872_p2601_identity_action_to_sigma_gamma_interface_audit.py"
OUT = ROOT / "generated" / "p2922_s1872_p2601_identity_action_to_sigma_gamma_interface_audit.json"
MD = ROOT / "generated" / "p2922_s1872_p2601_identity_action_to_sigma_gamma_interface_audit.md"


class P2922P2601IdentityActionToSigmaGammaInterfaceAuditTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        subprocess.run([sys.executable, str(SCRIPT)], cwd=ROOT, check=True, stdout=subprocess.PIPE, text=True)
        cls.payload = json.loads(OUT.read_text(encoding="utf-8"))

    def test_status_and_inputs(self):
        self.assertEqual(self.payload["status"], "P2922_P2601_IDENTITY_ACTION_TO_SIGMA_GAMMA_INTERFACE_AUDIT_NO_ACCEPTED_SOURCE")
        self.assertTrue(self.payload["acceptance_matrix"]["p2921_rechecked_verifier_exported"])
        self.assertTrue(self.payload["acceptance_matrix"]["p2601_identity_action_source_exported"])
        self.assertIsNotNone(self.payload["input_hashes"]["P2921"])
        self.assertIsNotNone(self.payload["input_hashes"]["P2601"])

    def test_interface_schema_and_source_facts(self):
        schema = self.payload["constructed_theoretical_objects"]["interface_schema"]
        facts = self.payload["constructed_theoretical_objects"]["p2601_source_facts"]
        self.assertEqual(schema["interface_name"], "P2601_Identity_Action_to_sigma_Gamma_Source_Interface")
        self.assertTrue(facts["hydrodynamic_identity_action_source_exported"])
        self.assertIn("y_1", facts["damping_amplitude_at_identity"])
        self.assertEqual(facts["remaining_keys_after_m2_and_unital"], ["prime_log_proportionality_source", "slope_value_or_prime_anchor_source"])

    def test_interface_candidate_counts(self):
        acc = self.payload["acceptance_matrix"]
        self.assertEqual(acc["interface_candidate_count"], 3)
        self.assertEqual(acc["accepted_interface_candidate_count"], 0)
        self.assertEqual(acc["nonzero_candidate_count"], 2)
        self.assertEqual(acc["gamma_lane_provenance_candidate_count"], 0)
        self.assertEqual(acc["iq_coupled_candidate_count"], 0)
        self.assertFalse(acc["strict_sigma_gamma_source_accepted"])

    def test_false_closure_exports_and_docs(self):
        flags = self.payload["decision"]["negative_export_flags"]
        for key in ["nonproxy_ltotal_exported", "eom_closure_exported", "hamiltonian_closure_exported", "bridge_closure_exported", "role_transfer_exported", "toe_closure_exported"]:
            self.assertFalse(flags[key])
        self.assertIn("P2922/S1872", MD.read_text(encoding="utf-8"))
        self.assertIn("P2922/S1872", (ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md").read_text(encoding="utf-8"))
        self.assertIn("P2922/S1872", (ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md").read_text(encoding="utf-8"))
        self.assertIn("P2922", (REPO / "AGENTS.md").read_text(encoding="utf-8"))


if __name__ == "__main__":
    unittest.main()
