import json, subprocess, sys, unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
SCRIPT = ROOT / "p3107_s2057_shannon_to_dimension_functor_source_law_audit.py"
OUT = ROOT / "generated" / "p3107_s2057_shannon_to_dimension_functor_source_law_audit.json"
MD = ROOT / "generated" / "p3107_s2057_shannon_to_dimension_functor_source_law_audit.md"

class P3107ShannonToDimensionFunctorSourceLawAuditTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        subprocess.run([sys.executable, str(SCRIPT)], cwd=ROOT, check=True, stdout=subprocess.PIPE, text=True)
        cls.payload = json.loads(OUT.read_text(encoding="utf-8"))

    def test_status_and_predecessor_hash(self):
        self.assertEqual(self.payload["status"], "P3107_SHANNON_TO_DIMENSION_FUNCTOR_SOURCE_LAW_OBSTRUCTION_BOUNDED_NO_GO")
        self.assertTrue(self.payload["input_hashes"]["P3106"])

    def test_certificate_counts(self):
        cert = self.payload["finite_certificate"]
        self.assertEqual(cert["p3106_accepted_internal_normalization_theorems"], 0)
        self.assertEqual(cert["information_category_rows"], 5)
        self.assertEqual(cert["alpha_geo_four_bit_cells"], 1)
        self.assertEqual(cert["additive_coherence_rows"], 15)
        self.assertEqual(cert["additive_coherence_failures"], 0)
        self.assertEqual(cert["dimension_functor_rows"], 35)
        self.assertEqual(cert["scale_invariance_rows"], 21)
        self.assertEqual(cert["functor_candidates"], 7)
        self.assertEqual(cert["required_gates"], 8)
        self.assertEqual(cert["candidate_gate_rows"], 56)
        self.assertEqual(cert["accepted_internal_shannon_to_dimension_source_laws"], 0)

    def test_functor_rows_show_exact_entropy_but_no_internal_units(self):
        objs = self.payload["constructed_theoretical_objects"]
        self.assertIn("self-coupled primordial information", objs["functor_audit_object"]["ontology"])
        self.assertTrue(any(row["is_alpha_geo_four_bit_cell"] for row in objs["information_category_rows"]))
        self.assertTrue(all(row["additive_shannon_coherence_passed"] for row in objs["additive_shannon_coherence_rows"]))
        self.assertTrue(any(row["nonzero_label"] for row in objs["dimension_functor_rows"]))
        self.assertTrue(all(not row["passes_without_import"] for row in objs["scale_invariance_rows"]))

    def test_negative_flags_and_recommendation(self):
        decision = self.payload["decision"]
        self.assertTrue(all(value is False for value in decision["negative_export_flags"].values()))
        self.assertIn("internal calibration-morphism candidate", decision["next_honest_step"])
        self.assertIn("P3107/S2057", MD.read_text(encoding="utf-8"))
        self.assertIn("P3107/S2057", (ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md").read_text(encoding="utf-8"))
        self.assertIn("P3107/S2057", (ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md").read_text(encoding="utf-8"))
        self.assertIn("P3107", (REPO / "AGENTS.md").read_text(encoding="utf-8"))

if __name__ == "__main__":
    unittest.main()
