import json, subprocess, sys, unittest
from pathlib import Path
ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
SCRIPT = ROOT / "p2973_s1923_incidence_unit_bearing_action_density_coupling_obstruction.py"
OUT = ROOT / "generated" / "p2973_s1923_incidence_unit_bearing_action_density_coupling_obstruction.json"
MD = ROOT / "generated" / "p2973_s1923_incidence_unit_bearing_action_density_coupling_obstruction.md"

class P2973IncidenceUnitBearingActionDensityCouplingObstructionTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        subprocess.run([sys.executable, str(SCRIPT)], cwd=ROOT, check=True, stdout=subprocess.PIPE, text=True)
        cls.payload = json.loads(OUT.read_text(encoding="utf-8"))

    def test_status_and_input(self):
        self.assertEqual(cls_status := self.payload["status"], "P2973_INCIDENCE_UNIT_BEARING_ACTION_DENSITY_COUPLING_OBSTRUCTION_NO_STRICT_EXPORT")
        self.assertIn("P2973", cls_status)
        self.assertIsNotNone(self.payload["input_hashes"]["P2972"])

    def test_coupling_certificate(self):
        cert = self.payload["coupling_certificate"]
        self.assertEqual(cert["slot_count"], 5)
        self.assertEqual(cert["incidence_weight_sum"], 9)
        self.assertEqual(cert["receiver_candidate_count"], 6)
        self.assertEqual(cert["available_receiver_candidate_count"], 5)
        self.assertEqual(cert["typed_reception_rows"], 5)
        self.assertEqual(cert["accepted_current_unit_bearing_couplings"], [])
        self.assertEqual(cert["acceptance_matrix_rows"], 64)
        self.assertEqual(cert["accepted_rows"], 1)

    def test_obligations_and_rows(self):
        obligations = {r["obligation"]: r["satisfied"] for r in self.payload["constructed_theoretical_objects"]["proof_obligation_rows"]}
        self.assertTrue(obligations["typed_incidence_reception_exists"])
        self.assertTrue(obligations["support_provenance_preserved_by_a_receiver"])
        self.assertFalse(obligations["nonconventional_unit_measure_source_exported"])
        self.assertFalse(obligations["strict_source_localizer_theorem_exported"])
        self.assertFalse(obligations["coupling_theorem_to_named_density_exported"])
        self.assertFalse(obligations["nonproxy_variational_readiness_exported"])
        self.assertFalse(obligations["accepted_current_unit_bearing_coupling"])
        rows = {r["candidate"]: r for r in self.payload["constructed_theoretical_objects"]["receiver_rows"]}
        self.assertTrue(rows["formal_slot_sum_density_L_incidence_formal"]["typed_incidence_reception"])
        self.assertFalse(rows["primitive_mean_9_5_density_import"]["typed_incidence_reception"])
        self.assertFalse(rows["completed_strict_unit_bearing_incidence_density_schema"]["accepted_current_unit_bearing_coupling"])

    def test_docs_and_nonpromotion(self):
        for value in self.payload["decision"]["negative_export_flags"].values():
            self.assertFalse(value)
        self.assertIn("P2973/S1923", MD.read_text(encoding="utf-8"))
        self.assertIn("P2973/S1923", (ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md").read_text(encoding="utf-8"))
        self.assertIn("P2973/S1923", (ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md").read_text(encoding="utf-8"))
        self.assertIn("P2973", (REPO / "AGENTS.md").read_text(encoding="utf-8"))

if __name__ == "__main__":
    unittest.main()
