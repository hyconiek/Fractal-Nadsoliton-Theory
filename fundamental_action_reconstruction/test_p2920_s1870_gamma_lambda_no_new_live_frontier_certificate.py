import json
import subprocess
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
SCRIPT = ROOT / "p2920_s1870_gamma_lambda_no_new_live_frontier_certificate.py"
OUT = ROOT / "generated" / "p2920_s1870_gamma_lambda_no_new_live_frontier_certificate.json"
MD = ROOT / "generated" / "p2920_s1870_gamma_lambda_no_new_live_frontier_certificate.md"


class P2920GammaLambdaNoNewLiveFrontierCertificateTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        subprocess.run([sys.executable, str(SCRIPT)], cwd=ROOT, check=True, stdout=subprocess.PIPE, text=True)
        cls.payload = json.loads(OUT.read_text(encoding="utf-8"))

    def test_status_and_p2919_input(self):
        self.assertEqual(self.payload["status"], "P2920_GAMMA_LAMBDA_NO_NEW_LIVE_FRONTIER_CERTIFICATE")
        self.assertTrue(self.payload["acceptance_matrix"]["p2919_rechecked_no_sigma_gamma_export"])
        self.assertIsNotNone(self.payload["input_hashes"]["P2919"])

    def test_state_map_counts(self):
        acc = self.payload["acceptance_matrix"]
        self.assertEqual(acc["state_map_row_count"], 9)
        self.assertEqual(acc["closed_lane_count"], 3)
        self.assertEqual(acc["finite_readiness_rows"], 9)
        self.assertEqual(acc["ltotal_unlocking_rows"], 0)

    def test_minimal_unlock_packet(self):
        packet = self.payload["constructed_theoretical_objects"]["minimal_unlock_packet"]
        self.assertEqual(packet["required_new_typed_object"], "Strict_sigma_Gamma_Action_Source_Map")
        self.assertIn("sigma_Gamma", packet["minimal_formula_shape"])
        self.assertEqual(len(packet["acceptance_tests"]), 5)
        self.assertEqual(self.payload["acceptance_matrix"]["minimal_unlock_object_count"], 1)

    def test_certificate_and_false_closure_flags(self):
        acc = self.payload["acceptance_matrix"]
        self.assertTrue(acc["no_new_live_frontier_certificate_exported"])
        self.assertFalse(acc["strict_sigma_gamma_source_exported"])
        self.assertFalse(acc["accepted_as_nonproxy_ltotal"])
        flags = self.payload["decision"]["negative_export_flags"]
        for key in ["nonproxy_ltotal_exported", "eom_closure_exported", "hamiltonian_closure_exported", "bridge_closure_exported", "role_transfer_exported", "toe_closure_exported"]:
            self.assertFalse(flags[key])

    def test_documents_updated(self):
        self.assertIn("P2920/S1870", MD.read_text(encoding="utf-8"))
        self.assertIn("P2920/S1870", (ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md").read_text(encoding="utf-8"))
        self.assertIn("P2920/S1870", (ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md").read_text(encoding="utf-8"))
        self.assertIn("P2920", (REPO / "AGENTS.md").read_text(encoding="utf-8"))


if __name__ == "__main__":
    unittest.main()
