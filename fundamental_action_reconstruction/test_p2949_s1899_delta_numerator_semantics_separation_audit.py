import json
import subprocess
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
SCRIPT = ROOT / "p2949_s1899_delta_numerator_semantics_separation_audit.py"
OUT = ROOT / "generated" / "p2949_s1899_delta_numerator_semantics_separation_audit.json"
MD = ROOT / "generated" / "p2949_s1899_delta_numerator_semantics_separation_audit.md"


class P2949DeltaNumeratorSemanticsSeparationAuditTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        subprocess.run([sys.executable, str(SCRIPT)], cwd=ROOT, check=True, stdout=subprocess.PIPE, text=True)
        cls.payload = json.loads(OUT.read_text(encoding="utf-8"))

    def test_status_and_inputs(self):
        self.assertEqual(self.payload["status"], "P2949_DELTA_NUMERATOR_SEMANTICS_SEPARATION_AUDIT_NO_STRICT_EXPORT")
        self.assertIsNotNone(self.payload["input_hashes"]["P2944"])
        self.assertIsNotNone(self.payload["input_hashes"]["P2948"])

    def test_delta_numerator_semantics_certificate(self):
        cert = self.payload["delta_numerator_semantics_certificate"]
        self.assertEqual(cert["identity_extension"], [1])
        self.assertEqual(cert["zero_extension"], [1])
        self.assertTrue(cert["identity_zero_extensions_equal"])
        self.assertTrue(cert["all_numerator_functionals_match_delta_4_5"])
        self.assertFalse(cert["strict_intensional_identity_deficit_semantics_exported"])
        self.assertFalse(cert["p2948_delta_numerator_premise_discharged"])

    def test_rows_acceptance_and_nonpromotion(self):
        obj = self.payload["constructed_theoretical_objects"]
        role_rows = obj["role_extension_rows"]
        self.assertEqual(sum(row["identity_role"] for row in role_rows), 1)
        self.assertEqual(sum(row["zero_role"] for row in role_rows), 1)
        self.assertTrue(all(row["roles_agree_on_node"] for row in role_rows))
        numerator_rows = obj["numerator_functional_rows"]
        self.assertEqual({row["delta"]["as_string"] for row in numerator_rows}, {"4/5"})
        self.assertTrue(all(not row["strict_semantics_selected"] for row in numerator_rows))
        acceptance = obj["acceptance_rows"]
        self.assertTrue(acceptance[0]["satisfied"])
        self.assertTrue(acceptance[1]["satisfied"])
        self.assertFalse(any(row["satisfied"] for row in acceptance[2:]))
        flags = self.payload["decision"]["negative_export_flags"]
        for key in ["strict_delta_numerator_semantics_exported", "strict_delta_eta_source_law_exported", "strict_beta_eta_coupling_theorem_exported", "strict_damping_beta_eta_source_packet_exported", "nonproxy_ltotal_exported", "eom_hamiltonian_exported", "bridge_closure_exported", "role_transfer_exported", "toe_closure_exported"]:
            self.assertFalse(flags[key])

    def test_docs_updated(self):
        self.assertIn("P2949/S1899", MD.read_text(encoding="utf-8"))
        self.assertIn("P2949/S1899", (ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md").read_text(encoding="utf-8"))
        self.assertIn("P2949/S1899", (ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md").read_text(encoding="utf-8"))
        self.assertIn("P2949", (REPO / "AGENTS.md").read_text(encoding="utf-8"))


if __name__ == "__main__":
    unittest.main()
