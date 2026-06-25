import json, subprocess, sys, unittest
from pathlib import Path
ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
SCRIPT = ROOT / "p3093_s2043_ward_identity_symmetry_current_effective_charge_obstruction_audit.py"
OUT = ROOT / "generated" / "p3093_s2043_ward_identity_symmetry_current_effective_charge_obstruction_audit.json"
MD = ROOT / "generated" / "p3093_s2043_ward_identity_symmetry_current_effective_charge_obstruction_audit.md"

class P3093WardIdentitySymmetryCurrentObstructionAuditTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        subprocess.run([sys.executable, str(SCRIPT)], cwd=ROOT, check=True, stdout=subprocess.PIPE, text=True)
        cls.payload = json.loads(OUT.read_text(encoding="utf-8"))

    def test_status_and_input_hash(self):
        self.assertEqual(self.payload["status"], "P3093_WARD_IDENTITY_SYMMETRY_CURRENT_EFFECTIVE_CHARGE_OBSTRUCTION_BOUNDED_NO_GO")
        self.assertIsNotNone(self.payload["input_hashes"]["P3092"])

    def test_finite_counts(self):
        cert = self.payload["finite_certificate"]
        self.assertEqual(cert["content_grep_lanes"], 3)
        self.assertGreater(cert["content_grep_hits"], 0)
        self.assertEqual(cert["p3092_accepted_nonimported_renormalization_scale_flow_sources"], 0)
        self.assertEqual(cert["symmetry_commutator_rows"], 24)
        self.assertEqual(cert["symmetry_commutator_failures"], 0)
        self.assertEqual(cert["symmetry_rows_with_current_or_charge_semantics"], 0)
        self.assertEqual(cert["source_variation_rows"], 12)
        self.assertEqual(cert["source_variation_balance_failures"], 0)
        self.assertEqual(cert["source_variation_rows_with_ward_identity"], 0)
        self.assertEqual(cert["spectral_charge_label_rows"], 12)
        self.assertEqual(cert["spectral_charge_rows_with_physical_normalization"], 0)
        self.assertEqual(cert["spectral_charge_rows_with_empirical_readout"], 0)
        self.assertEqual(cert["formal_discrete_current_rows"], 12)
        self.assertEqual(cert["formal_current_rows_with_physical_current_law"], 0)
        self.assertEqual(cert["charge_candidates"], 5)
        self.assertEqual(cert["candidate_gate_rows"], 40)
        self.assertEqual(cert["accepted_nonimported_ward_current_effective_charge_sources"], 0)
        self.assertEqual(cert["satisfied_proof_obligations"], 4)

    def test_witnesses_and_candidates(self):
        objs = self.payload["constructed_theoretical_objects"]
        self.assertTrue(all(row["laplacian_commutator_identity_witness"] for row in objs["symmetry_commutator_rows"]))
        self.assertTrue(all(row["finite_source_variation_balance_witness"] for row in objs["source_variation_balance_rows"]))
        self.assertEqual(len(objs["spectral_charge_label_rows"]), 12)
        self.assertTrue(all(row["current_is_formal_profile_gradient"] for row in objs["formal_discrete_current_rows"]))
        aggs = objs["candidate_aggregate_certificate"]
        self.assertTrue(all(not row["accepted_nonimported_ward_current_effective_charge_source"] for row in aggs))
        by_id = {row["candidate"]: row for row in aggs}
        self.assertEqual(by_id["z12_graph_automorphism_commutant"]["passed_gates"], 2)
        self.assertEqual(by_id["finite_source_variation_balance"]["passed_gates"], 3)
        self.assertEqual(by_id["imported_empirical_electric_charge_template"]["passed_gates"], 7)

    def test_flags_recommendation_and_docs(self):
        for value in self.payload["decision"]["negative_export_flags"].values():
            self.assertFalse(value)
        self.assertTrue(self.payload["decision"]["positive_scoped_flags"]["finite_graph_symmetry_commutators_computed"])
        self.assertIn("stress-energy/metric-response", self.payload["decision"]["next_honest_step"])
        self.assertIn("P3093/S2043", MD.read_text(encoding="utf-8"))
        self.assertIn("P3093/S2043", (ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md").read_text(encoding="utf-8"))
        self.assertIn("P3093/S2043", (ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md").read_text(encoding="utf-8"))
        self.assertIn("P3093", (REPO / "AGENTS.md").read_text(encoding="utf-8"))

if __name__ == "__main__":
    unittest.main()
