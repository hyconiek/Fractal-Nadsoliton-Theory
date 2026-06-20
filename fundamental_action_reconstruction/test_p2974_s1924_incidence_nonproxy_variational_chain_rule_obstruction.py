import json, subprocess, sys, unittest
from pathlib import Path
ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
SCRIPT = ROOT / "p2974_s1924_incidence_nonproxy_variational_chain_rule_obstruction.py"
OUT = ROOT / "generated" / "p2974_s1924_incidence_nonproxy_variational_chain_rule_obstruction.json"
MD = ROOT / "generated" / "p2974_s1924_incidence_nonproxy_variational_chain_rule_obstruction.md"

class P2974IncidenceNonproxyVariationalChainRuleObstructionTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        subprocess.run([sys.executable, str(SCRIPT)], cwd=ROOT, check=True, stdout=subprocess.PIPE, text=True)
        cls.payload = json.loads(OUT.read_text(encoding="utf-8"))

    def test_status_and_input(self):
        self.assertEqual(self.payload["status"], "P2974_INCIDENCE_NONPROXY_VARIATIONAL_CHAIN_RULE_OBSTRUCTION_NO_STRICT_EXPORT")
        self.assertIsNotNone(self.payload["input_hashes"]["P2973"])

    def test_chain_rule_certificate(self):
        cert = self.payload["chain_rule_certificate"]
        self.assertEqual(cert["slot_count"], 5)
        self.assertEqual(cert["weight_sum"], 9)
        self.assertEqual(cert["nonzero_hessian_entries"], 5)
        self.assertEqual(cert["off_diagonal_entries"], 0)
        self.assertEqual(cert["rank_over_Q"], 5)
        self.assertEqual(cert["candidate_count"], 5)
        self.assertEqual(cert["accepted_current_nonproxy_chain_rules"], [])
        self.assertEqual(cert["acceptance_matrix_rows"], 128)
        self.assertEqual(cert["accepted_rows"], 1)

    def test_obligations_and_rows(self):
        obligations = {r["obligation"]: r["satisfied"] for r in self.payload["constructed_theoretical_objects"]["proof_obligation_rows"]}
        self.assertTrue(obligations["finite_jacobian_constructed"])
        self.assertTrue(obligations["incidence_slots_preserved"])
        self.assertFalse(obligations["strict_field_variable_provenance_exported"])
        self.assertFalse(obligations["unit_measure_coupled_density_exported"])
        self.assertFalse(obligations["source_localizer_available"])
        self.assertFalse(obligations["boundary_integration_theorem_exported"])
        self.assertFalse(obligations["continuum_nonproxy_lift_exported"])
        self.assertFalse(obligations["accepted_current_nonproxy_chain_rule"])
        jac = self.payload["constructed_theoretical_objects"]["formal_quadratic_jacobian"]
        self.assertEqual(jac["weights"], [1, 2, 2, 2, 2])
        self.assertEqual(jac["hessian_matrix"][0], [1, 0, 0, 0, 0])
        rows = {r["candidate"]: r for r in self.payload["constructed_theoretical_objects"]["chain_rule_candidate_rows"]}
        self.assertFalse(rows["completed_strict_incidence_nonproxy_chain_rule_schema"]["accepted_current_nonproxy_chain_rule"])

    def test_docs_and_nonpromotion(self):
        for value in self.payload["decision"]["negative_export_flags"].values():
            self.assertFalse(value)
        self.assertIn("P2974/S1924", MD.read_text(encoding="utf-8"))
        self.assertIn("P2974/S1924", (ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md").read_text(encoding="utf-8"))
        self.assertIn("P2974/S1924", (ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md").read_text(encoding="utf-8"))
        self.assertIn("P2974", (REPO / "AGENTS.md").read_text(encoding="utf-8"))

if __name__ == "__main__":
    unittest.main()
