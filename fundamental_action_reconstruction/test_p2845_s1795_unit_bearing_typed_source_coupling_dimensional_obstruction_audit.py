import json
import subprocess
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
SCRIPT = ROOT / "p2845_s1795_unit_bearing_typed_source_coupling_dimensional_obstruction_audit.py"
JSON_PATH = ROOT / "generated" / "p2845_s1795_unit_bearing_typed_source_coupling_dimensional_obstruction_audit.json"
MD_PATH = ROOT / "generated" / "p2845_s1795_unit_bearing_typed_source_coupling_dimensional_obstruction_audit.md"


class P2845UnitBearingTypedSourceCouplingDimensionalAuditTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        subprocess.run([sys.executable, str(SCRIPT)], check=True, cwd=ROOT)
        cls.payload = json.loads(JSON_PATH.read_text(encoding="utf-8"))

    def test_status_and_input(self):
        self.assertEqual(
            self.payload["status"],
            "P2845_UNIT_BEARING_TYPED_SOURCE_COUPLING_DIMENSIONAL_OBSTRUCTION_NO_CLOSURE",
        )
        audit = self.payload["unit_bearing_typed_source_coupling_audit"]
        self.assertEqual(
            audit["input_statuses_rechecked"]["P2844"],
            "P2844_CLOSURE_GATE_PRIME_IMPLICANT_OBLIGATION_MATRIX_NO_CLOSURE",
        )

    def test_dimension_conventions_and_rows(self):
        audit = self.payload["unit_bearing_typed_source_coupling_audit"]
        self.assertEqual(audit["dimension_conventions"]["spacetime_dimension"], 4)
        self.assertEqual(audit["dimension_conventions"]["lagrangian_density_mass_dimension"], 4)
        self.assertEqual(len(audit["candidate_rows"]), 5)
        required_dims = {row["candidate"]: row["required_coupling_dimension"]["display"] for row in audit["candidate_rows"]}
        self.assertEqual(required_dims["dimensionless_graph_times_dimension4_scalar_density"], "0")
        self.assertEqual(required_dims["graph_mass_term_scalar_square"], "2")
        self.assertEqual(required_dims["strict_kernel_modulated_density"], "0")

    def test_all_candidates_formally_balanced_but_rejected(self):
        audit = self.payload["unit_bearing_typed_source_coupling_audit"]
        self.assertEqual(audit["formal_dimension_balance_count"], 5)
        self.assertEqual(audit["target_independent_units_count"], 0)
        self.assertEqual(audit["accepted_candidate_count"], 0)
        for row in audit["candidate_rows"]:
            self.assertTrue(row["premises"]["dimension_balanced_units"])
            self.assertFalse(row["premises"]["target_independent_units"])
            self.assertFalse(row["accepted_as_unit_bearing_typed_source_coupling"])

    def test_blockers_and_hamiltonian_boundary(self):
        audit = self.payload["unit_bearing_typed_source_coupling_audit"]
        self.assertEqual(audit["blocker_histogram"]["localization_pullback"], 5)
        self.assertEqual(audit["blocker_histogram"]["coupling_coefficient_rule"], 5)
        h_row = next(row for row in audit["candidate_rows"] if row["candidate"] == "hamiltonian_potential_placeholder")
        self.assertIn("typed_target_codomain", h_row["missing_premises"])
        self.assertIn("nonproxy_ltotal_term", h_row["missing_premises"])

    def test_acceptance_and_negative_exports(self):
        acceptance = self.payload["acceptance_matrix"]
        self.assertTrue(acceptance["accepted_as_unit_bearing_coupling_obstruction_audit"])
        self.assertFalse(acceptance["exports_unit_bearing_typed_source_coupling"])
        flags = self.payload["decision"]["negative_export_flags"]
        self.assertFalse(flags["unit_bearing_typed_source_coupling_exported"])
        self.assertFalse(flags["nonproxy_ltotal_term_exported"])
        self.assertFalse(flags["eom_closure_exported"])
        self.assertFalse(flags["hamiltonian_closure_exported"])
        self.assertFalse(flags["toe_closure_exported"])

    def test_documents_updated(self):
        self.assertIn("P2845/S1795", MD_PATH.read_text(encoding="utf-8"))
        self.assertIn("P2845/S1795", (ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md").read_text(encoding="utf-8"))
        self.assertIn("P2845/S1795", (ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md").read_text(encoding="utf-8"))
        self.assertIn("P2845", (REPO / "AGENTS.md").read_text(encoding="utf-8"))


if __name__ == "__main__":
    unittest.main()
