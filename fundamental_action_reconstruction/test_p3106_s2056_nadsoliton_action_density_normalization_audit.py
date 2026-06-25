import json, subprocess, sys, unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
SCRIPT = ROOT / "p3106_s2056_nadsoliton_action_density_normalization_audit.py"
OUT = ROOT / "generated" / "p3106_s2056_nadsoliton_action_density_normalization_audit.json"
MD = ROOT / "generated" / "p3106_s2056_nadsoliton_action_density_normalization_audit.md"

class P3106NadsolitonActionDensityNormalizationAuditTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        subprocess.run([sys.executable, str(SCRIPT)], cwd=ROOT, check=True, stdout=subprocess.PIPE, text=True)
        cls.payload = json.loads(OUT.read_text(encoding="utf-8"))

    def test_status_and_predecessor_hash(self):
        self.assertEqual(self.payload["status"], "P3106_NADSOLITON_ACTION_DENSITY_NORMALIZATION_OBSTRUCTION_BOUNDED_NO_GO")
        self.assertTrue(self.payload["input_hashes"]["P3105"])

    def test_certificate_counts(self):
        cert = self.payload["finite_certificate"]
        self.assertEqual(cert["p3105_accepted_nonimported_unit_sources"], 0)
        self.assertEqual(cert["shannon_cell_rows"], 4)
        self.assertEqual(cert["four_bit_alpha_cells"], 1)
        self.assertEqual(cert["bit_to_action_dimension_laws"], 0)
        self.assertEqual(cert["action_density_rows"], 16)
        self.assertEqual(cert["density_rows_with_dimension_assignment"], 0)
        self.assertEqual(cert["scale_quotient_rows"], 4)
        self.assertEqual(cert["nonconventional_scale_representatives"], 0)
        self.assertEqual(cert["normalization_candidates"], 6)
        self.assertEqual(cert["required_gates"], 7)
        self.assertEqual(cert["candidate_gate_rows"], 42)
        self.assertEqual(cert["accepted_internal_normalization_theorems"], 0)

    def test_constructed_rows_are_internal_but_unitless(self):
        objs = self.payload["constructed_theoretical_objects"]
        self.assertIn("self-coupled primordial information", objs["normalization_audit_object"]["ontology"])
        self.assertTrue(any(row["is_four_bit_alpha_cell"] for row in objs["shannon_four_bit_cell_rows"]))
        self.assertTrue(all(not row["bit_to_action_dimension_law_exported"] for row in objs["shannon_four_bit_cell_rows"]))
        self.assertTrue(all(row["self_coupled_density_defined"] for row in objs["action_density_rows"]))
        self.assertTrue(all(not row["dimension_assignment_exported"] for row in objs["action_density_rows"]))
        self.assertTrue(all(row["positive_scale_quotient_constructed"] for row in objs["scale_quotient_rows"]))
        self.assertTrue(all(not row["canonical_representative_selected_nonconventionally"] for row in objs["scale_quotient_rows"]))

    def test_negative_flags_and_recommendation(self):
        decision = self.payload["decision"]
        self.assertTrue(all(value is False for value in decision["negative_export_flags"].values()))
        self.assertIn("Shannon-to-dimension functor", decision["next_honest_step"])
        self.assertIn("P3106/S2056", MD.read_text(encoding="utf-8"))
        self.assertIn("P3106/S2056", (ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md").read_text(encoding="utf-8"))
        self.assertIn("P3106/S2056", (ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md").read_text(encoding="utf-8"))
        self.assertIn("P3106", (REPO / "AGENTS.md").read_text(encoding="utf-8"))

if __name__ == "__main__":
    unittest.main()
