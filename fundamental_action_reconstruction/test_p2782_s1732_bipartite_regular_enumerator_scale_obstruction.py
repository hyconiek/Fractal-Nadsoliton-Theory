import json
import subprocess
import unittest
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]
ROOT = REPO / "fundamental_action_reconstruction"
SCRIPT = ROOT / "p2782_s1732_bipartite_regular_enumerator_scale_obstruction.py"
OUT = ROOT / "generated" / "p2782_s1732_bipartite_regular_enumerator_scale_obstruction.json"
MD = ROOT / "generated" / "p2782_s1732_bipartite_regular_enumerator_scale_obstruction.md"


class P2782BipartiteRegularEnumeratorScaleObstructionTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        subprocess.run(["python", str(SCRIPT)], cwd=ROOT, check=True)
        cls.payload = json.loads(OUT.read_text(encoding="utf-8"))

    def test_status_and_input(self):
        self.assertEqual(self.payload["status"], "P2782_BIPARTITE_REGULAR_ENUMERATOR_SCALE_OBSTRUCTION_NO_CLOSURE")
        self.assertEqual(self.payload["input_statuses"]["P2781"], "P2781_ENUMERATED_TWO_LAYER_C8_SPECTRUM_COLLISION_AUDIT_NO_CLOSURE")
        self.assertIn("naive canonical enumerator", self.payload["audited_question"])

    def test_exact_dp_count_blocks_naive_enumerator(self):
        witness = self.payload["bipartite_regular_enumerator_scale_witness"]
        count_row = witness["count_row"]
        self.assertEqual(count_row["left_size"], 8)
        self.assertEqual(count_row["right_size"], 8)
        self.assertEqual(count_row["degree"], 4)
        self.assertEqual(count_row["row_mask_count"], 70)
        self.assertEqual(count_row["labeled_bipartite_regular_matrix_count"], 116_963_796_250)
        self.assertTrue(count_row["matches_expected_regression_value"])
        self.assertEqual(witness["row_column_relabeling_group_size"], 1_625_702_400)
        self.assertTrue(witness["naive_enumerator_blocked_in_repo"])
        self.assertTrue(witness["connectivity_filter_not_yet_applied"])
        self.assertTrue(witness["graph_isomorphism_quotient_not_yet_applied"])

    def test_acceptance_is_obstruction_not_closure(self):
        acceptance = self.payload["acceptance_matrix"]
        self.assertTrue(acceptance["facts"]["exact_dp_count_performed"])
        self.assertTrue(acceptance["facts"]["fixed_bipartition_subproblem_exceeds_naive_in_repo_scale"])
        self.assertTrue(acceptance["accepted_as_enumerator_scale_obstruction"])
        self.assertFalse(acceptance["accepted_as_canonical_16node_4regular_enumerator"])
        self.assertFalse(acceptance["accepted_as_strict_spectral_source_law"])
        self.assertFalse(acceptance["accepted_as_canonical_nadsoliton_geometry_source"])
        self.assertIn("canonical_generation_algorithm_or_external_certificate_supplied", acceptance["missing_criteria"])
        self.assertTrue(all(value is False for value in self.payload["decision"]["negative_export_flags"].values()))

    def test_documentation_and_recommendation(self):
        self.assertIn("canonical-generation theorem/tool certificate", self.payload["decision"]["next_honest_step"])
        self.assertIn("P2697-P2782", self.payload["decision"]["next_honest_step"])
        self.assertIn("P2782/S1732", MD.read_text(encoding="utf-8"))
        self.assertIn("P2782/S1732", (ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md").read_text(encoding="utf-8"))
        self.assertIn("P2782/S1732", (ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md").read_text(encoding="utf-8"))
        self.assertIn("P2782", (REPO / "AGENTS.md").read_text(encoding="utf-8"))


if __name__ == "__main__":
    unittest.main()
