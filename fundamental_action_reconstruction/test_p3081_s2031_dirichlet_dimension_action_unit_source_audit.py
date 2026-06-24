import json, subprocess, sys, unittest
from pathlib import Path
ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
SCRIPT = ROOT / "p3081_s2031_dirichlet_dimension_action_unit_source_audit.py"
OUT = ROOT / "generated" / "p3081_s2031_dirichlet_dimension_action_unit_source_audit.json"
MD = ROOT / "generated" / "p3081_s2031_dirichlet_dimension_action_unit_source_audit.md"

class P3081DirichletDimensionActionUnitSourceAuditTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        subprocess.run([sys.executable, str(SCRIPT)], cwd=ROOT, check=True, stdout=subprocess.PIPE, text=True)
        cls.payload = json.loads(OUT.read_text(encoding="utf-8"))

    def test_status_and_input_hash(self):
        self.assertEqual(self.payload["status"], "P3081_DIRICHLET_DIMENSION_ACTION_UNIT_SOURCE_BOUNDED_NO_GO")
        self.assertIsNotNone(self.payload["input_hashes"]["P3080"])

    def test_finite_counts(self):
        cert = self.payload["finite_certificate"]
        self.assertEqual(cert["content_grep_lanes"], 3)
        self.assertGreater(cert["content_grep_hits"], 0)
        self.assertEqual(cert["p3080_accepted_standard_physics_interface_objects"], 0)
        self.assertEqual(cert["binary_profile_rows"], 4096)
        self.assertEqual(cert["constant_profile_rows"], 2)
        self.assertEqual(cert["nonconstant_profile_rows"], 4094)
        self.assertEqual(cert["energy_spectrum_rows"], 7)
        self.assertEqual(cert["minimum_nonzero_dirichlet_energy"], 1.0)
        self.assertEqual(cert["maximum_dirichlet_energy"], 6.0)
        self.assertEqual(cert["unit_source_candidates"], 6)
        self.assertEqual(cert["candidate_gate_rows"], 36)
        self.assertEqual(cert["accepted_nonimported_dimension_action_unit_sources"], 0)
        self.assertEqual(cert["satisfied_proof_obligations"], 4)

    def test_spectrum_and_candidate_aggregates(self):
        spectrum = self.payload["constructed_theoretical_objects"]["binary_dirichlet_energy_spectrum"]
        self.assertEqual([row["boundary_edges"] for row in spectrum], [0, 2, 4, 6, 8, 10, 12])
        self.assertEqual(sum(row["profile_count"] for row in spectrum), 4096)
        aggs = self.payload["constructed_theoretical_objects"]["candidate_aggregate_certificate"]
        self.assertTrue(all(not row["accepted_nonimported_dimension_action_unit_source"] for row in aggs))
        by_id = {row["candidate"]: row for row in aggs}
        self.assertEqual(by_id["cycle_cardinality_twelve"]["passed_gates"], 4)
        self.assertEqual(by_id["imported_hbar_c_lattice_spacing_template"]["passed_gates"], 5)

    def test_flags_recommendation_and_docs(self):
        for value in self.payload["decision"]["negative_export_flags"].values():
            self.assertFalse(value)
        self.assertTrue(self.payload["decision"]["positive_scoped_flags"]["binary_dirichlet_spectrum_computed"])
        self.assertIn("continuum-limit functor", self.payload["decision"]["next_honest_step"])
        self.assertIn("P3081/S2031", MD.read_text(encoding="utf-8"))
        self.assertIn("P3081/S2031", (ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md").read_text(encoding="utf-8"))
        self.assertIn("P3081/S2031", (ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md").read_text(encoding="utf-8"))
        self.assertIn("P3081", (REPO / "AGENTS.md").read_text(encoding="utf-8"))

if __name__ == "__main__":
    unittest.main()
