import json
import subprocess
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
SCRIPT = ROOT / "p2932_s1882_aut_z12_invariant_additive_value_source_no_go.py"
OUT = ROOT / "generated" / "p2932_s1882_aut_z12_invariant_additive_value_source_no_go.json"
MD = ROOT / "generated" / "p2932_s1882_aut_z12_invariant_additive_value_source_no_go.md"


class P2932AutZ12InvariantAdditiveValueSourceNoGoTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        subprocess.run([sys.executable, str(SCRIPT)], cwd=ROOT, check=True, stdout=subprocess.PIPE, text=True)
        cls.payload = json.loads(OUT.read_text(encoding="utf-8"))

    def test_status_and_input(self):
        self.assertEqual(self.payload["status"], "P2932_AUT_Z12_INVARIANT_ADDITIVE_VALUE_SOURCE_NO_GO_ZERO_ONLY")
        self.assertIsNotNone(self.payload["input_hashes"]["P2931"])

    def test_exact_rref_no_go_certificate(self):
        cert = self.payload["no_go_certificate"]
        self.assertEqual(cert["equation_count"], 74)
        self.assertEqual(cert["product_equation_count"], 29)
        self.assertEqual(cert["aut_invariance_equation_count"], 44)
        self.assertEqual(cert["rank"], 11)
        self.assertEqual(cert["nullity"], 0)
        self.assertTrue(cert["zero_function_only"])
        self.assertFalse(cert["nonzero_strict_prime_log_value_source_exported"])

    def test_candidate_space_and_no_closure_flags(self):
        candidate = self.payload["constructed_theoretical_objects"]["candidate_space"]
        self.assertEqual(candidate["name"], "AutZ12_Invariant_Additive_Prime_Value_Source_Space")
        self.assertEqual(len(candidate["obligations"]), 3)
        flags = self.payload["decision"]["negative_export_flags"]
        for key in ["strict_prime_log_value_source_exported", "strict_delta_eta_source_exported", "strict_beta_eta_coupling_theorem_exported", "strict_damping_beta_eta_source_packet_exported", "nonproxy_ltotal_exported", "eom_hamiltonian_exported", "bridge_closure_exported", "role_transfer_exported", "toe_closure_exported"]:
            self.assertFalse(flags[key])

    def test_docs_updated(self):
        self.assertIn("P2932/S1882", MD.read_text(encoding="utf-8"))
        self.assertIn("P2932/S1882", (ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md").read_text(encoding="utf-8"))
        self.assertIn("P2932/S1882", (ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md").read_text(encoding="utf-8"))
        self.assertIn("P2932", (REPO / "AGENTS.md").read_text(encoding="utf-8"))


if __name__ == "__main__":
    unittest.main()
