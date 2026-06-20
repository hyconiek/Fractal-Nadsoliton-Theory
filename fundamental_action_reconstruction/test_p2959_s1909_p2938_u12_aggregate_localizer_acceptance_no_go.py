import json
import subprocess
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
SCRIPT = ROOT / "p2959_s1909_p2938_u12_aggregate_localizer_acceptance_no_go.py"
OUT = ROOT / "generated" / "p2959_s1909_p2938_u12_aggregate_localizer_acceptance_no_go.json"
MD = ROOT / "generated" / "p2959_s1909_p2938_u12_aggregate_localizer_acceptance_no_go.md"


class P2959P2938U12AggregateLocalizerAcceptanceNoGoTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        subprocess.run([sys.executable, str(SCRIPT)], cwd=ROOT, check=True, stdout=subprocess.PIPE, text=True)
        cls.payload = json.loads(OUT.read_text(encoding="utf-8"))

    def test_status_and_inputs(self):
        self.assertEqual(self.payload["status"], "P2959_P2938_U12_AGGREGATE_LOCALIZER_ACCEPTANCE_NO_GO_NO_STRICT_EXPORT")
        self.assertIsNotNone(self.payload["input_hashes"]["P2953"])
        self.assertIsNotNone(self.payload["input_hashes"]["P2958"])

    def test_candidate_lattice_and_predicates(self):
        cert = self.payload["localizer_certificate"]
        self.assertEqual(cert["candidate_lattice_row_count"], 15)
        self.assertEqual(cert["exact_target_pair"], [1, 1])
        self.assertEqual(cert["exact_target_vector"], [1, 2, 2, 2, 2])
        self.assertNotIn("sum_9", cert["predicates_selecting_exact_target_only"])
        self.assertIn("primitive_equal_summand", cert["predicates_selecting_exact_target_only"])
        self.assertFalse(cert["all_exact_selecting_predicates_strict_exported"])
        predicates = {row["predicate"]: row for row in self.payload["constructed_theoretical_objects"]["localizer_predicate_rows"]}
        self.assertGreater(predicates["positive_coordinates"]["selected_count"], 1)
        self.assertFalse(predicates["sum_9"]["selects_exact_target_pair_only"])
        self.assertEqual(predicates["sum_9"]["selected_count"], 2)
        self.assertFalse(predicates["sum_9"]["strict_nadsoliton_localizer_exported"])

    def test_acceptance_obligations_and_matrix(self):
        obligations = {row["obligation"]: row for row in self.payload["constructed_theoretical_objects"]["acceptance_obligation_rows"]}
        self.assertTrue(obligations["finite_localizer_lattice_constructed"]["satisfied"])
        self.assertTrue(obligations["target_pair_selectable_by_some_predicate"]["satisfied"])
        self.assertFalse(obligations["predicate_not_target_coded_or_conventional"]["satisfied"])
        self.assertFalse(obligations["strict_nadsoliton_functor_exports_predicate"]["satisfied"])
        self.assertFalse(obligations["strict_equal_weight_source_theorem_exported"]["satisfied"])
        self.assertFalse(obligations["downstream_beta_unit_coupling_exported"]["satisfied"])
        matrix = self.payload["constructed_theoretical_objects"]["finite_acceptance_matrix"]
        self.assertEqual(len(matrix), 64)
        self.assertEqual(sum(1 for row in matrix if row["accepts_strict_u12_aggregate_localizer"]), 1)

    def test_nonpromotion_and_docs(self):
        cert = self.payload["localizer_certificate"]
        self.assertFalse(cert["p2958_functor_localizer_obligation_discharged"])
        for value in self.payload["decision"]["negative_export_flags"].values():
            self.assertFalse(value)
        self.assertIn("P2959/S1909", MD.read_text(encoding="utf-8"))
        self.assertIn("P2959/S1909", (ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md").read_text(encoding="utf-8"))
        self.assertIn("P2959/S1909", (ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md").read_text(encoding="utf-8"))
        self.assertIn("P2959", (REPO / "AGENTS.md").read_text(encoding="utf-8"))


if __name__ == "__main__":
    unittest.main()
