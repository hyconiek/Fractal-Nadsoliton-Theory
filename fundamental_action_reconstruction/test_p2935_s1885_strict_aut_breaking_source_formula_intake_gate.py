import json
import subprocess
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
SCRIPT = ROOT / "p2935_s1885_strict_aut_breaking_source_formula_intake_gate.py"
OUT = ROOT / "generated" / "p2935_s1885_strict_aut_breaking_source_formula_intake_gate.json"
MD = ROOT / "generated" / "p2935_s1885_strict_aut_breaking_source_formula_intake_gate.md"


class P2935StrictAutBreakingSourceFormulaIntakeGateTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        subprocess.run([sys.executable, str(SCRIPT)], cwd=ROOT, check=True, stdout=subprocess.PIPE, text=True)
        cls.payload = json.loads(OUT.read_text(encoding="utf-8"))

    def test_status_and_input(self):
        self.assertEqual(self.payload["status"], "P2935_STRICT_AUT_BREAKING_SOURCE_FORMULA_INTAKE_GATE_NO_ACCEPTED_FORMULA")
        self.assertIsNotNone(self.payload["input_hashes"]["P2934"])

    def test_intake_certificate(self):
        cert = self.payload["intake_certificate"]
        self.assertEqual(cert["obligation_count"], 5)
        self.assertEqual(cert["status_table_row_count"], 32)
        self.assertEqual(cert["accepting_row_count"], 1)
        self.assertEqual(cert["current_algebraic_readiness_only_row_count"], 1)
        self.assertEqual(cert["candidate_formula_class_count"], 5)
        self.assertEqual(cert["accepted_current_candidate_count"], 0)
        self.assertTrue(cert["no_new_live_frontier_certificate_exported"])

    def test_gate_and_no_closure_flags(self):
        gate = self.payload["constructed_theoretical_objects"]["intake_gate"]
        self.assertEqual(gate["name"], "Strict_AutBreaking_PrimeCoordinate_Source_Formula_Intake_Gate")
        self.assertEqual(len(gate["obligations"]), 5)
        flags = self.payload["decision"]["negative_export_flags"]
        for key in ["strict_aut_breaking_prime_coordinate_source_law_exported", "strict_prime_log_value_source_exported", "strict_delta_eta_source_exported", "strict_beta_eta_coupling_theorem_exported", "strict_damping_beta_eta_source_packet_exported", "nonproxy_ltotal_exported", "eom_hamiltonian_exported", "bridge_closure_exported", "role_transfer_exported", "toe_closure_exported"]:
            self.assertFalse(flags[key])

    def test_docs_updated(self):
        self.assertIn("P2935/S1885", MD.read_text(encoding="utf-8"))
        self.assertIn("P2935/S1885", (ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md").read_text(encoding="utf-8"))
        self.assertIn("P2935/S1885", (ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md").read_text(encoding="utf-8"))
        self.assertIn("P2935", (REPO / "AGENTS.md").read_text(encoding="utf-8"))


if __name__ == "__main__":
    unittest.main()
