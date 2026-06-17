import json
import subprocess
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
SCRIPT = ROOT / "p2831_s1781_p2815_digest_source_law_coupling_obligation_matrix.py"
JSON_PATH = ROOT / "generated" / "p2831_s1781_p2815_digest_source_law_coupling_obligation_matrix.json"
MD_PATH = ROOT / "generated" / "p2831_s1781_p2815_digest_source_law_coupling_obligation_matrix.md"


class P2831P2815DigestSourceLawCouplingObligationMatrixTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        subprocess.run([sys.executable, str(SCRIPT)], check=True, cwd=ROOT)
        cls.payload = json.loads(JSON_PATH.read_text(encoding="utf-8"))

    def test_status_and_inputs(self):
        self.assertEqual(self.payload["status"], "P2831_P2815_DIGEST_SOURCE_LAW_COUPLING_OBLIGATION_NO_GO_NO_CLOSURE")
        self.assertEqual(self.payload["input_statuses"]["P2818"], "P2818_LOCAL_EDGE_VARIATIONAL_RESPONSE_ENERGY_OBSTRUCTED_NO_COUPLING_NO_CLOSURE")
        self.assertEqual(self.payload["input_statuses"]["P2830"], "P2830_P2815_DIGEST_FULL_16828_CARRIER_SEPARATION_WITNESS_NO_SOURCE_LAW_NO_CLOSURE")

    def test_obligation_counts(self):
        audit = self.payload["p2815_digest_source_law_coupling_obligation_matrix"]
        self.assertEqual(audit["decoded_graph_count"], 16828)
        self.assertEqual(audit["p2830_full_carrier_rows"], 16828)
        self.assertEqual(audit["p2830_full_carrier_digest_classes"], 16828)
        self.assertEqual(audit["p2830_full_carrier_digest_collisions"], 0)
        self.assertEqual(audit["graph_index_min"], 0)
        self.assertEqual(audit["graph_index_max"], 16827)
        self.assertTrue(audit["graph_index_coverage_is_exact"])
        self.assertEqual(audit["accepted_candidate_count"], 0)
        self.assertEqual(audit["rejected_candidate_count"], 4)

    def test_candidate_rows(self):
        rows = {row["candidate"]: row for row in self.payload["p2815_digest_source_law_coupling_obligation_matrix"]["candidate_obligation_rows"]}
        self.assertEqual(set(rows), {"Q_digest_sha256_label", "Q_digest_lex_rank", "Q_digest_prefix_16", "Q_digest_int_mod_4096"})
        for name in ["Q_digest_sha256_label", "Q_digest_lex_rank", "Q_digest_prefix_16"]:
            self.assertEqual(rows[name]["value_count"], 16828)
            self.assertEqual(rows[name]["unique_value_count"], 16828)
            self.assertEqual(rows[name]["collision_class_count"], 0)
            self.assertFalse(rows[name]["accepted_as_source_law_coupling"])
            self.assertIn("non_label_graph_formula_exported", rows[name]["missing_for_promotion"])
            self.assertIn("target_independent_units_or_normalization_exported", rows[name]["missing_for_promotion"])
            self.assertIn("variational_derivative_exported", rows[name]["missing_for_promotion"])
            self.assertIn("typed_graph_to_K_or_Ltotal_coupling_theorem_exported", rows[name]["missing_for_promotion"])
        self.assertEqual(rows["Q_digest_int_mod_4096"]["unique_value_count"], 4014)
        self.assertEqual(rows["Q_digest_int_mod_4096"]["collision_class_count"], 3729)
        self.assertEqual(rows["Q_digest_int_mod_4096"]["max_class_size"], 13)
        self.assertIn("separates_full_carrier", rows["Q_digest_int_mod_4096"]["missing_for_promotion"])

    def test_acceptance_boundaries(self):
        matrix = self.payload["acceptance_matrix"]
        self.assertTrue(matrix["facts"]["p2818_local_response_rejected"])
        self.assertTrue(matrix["facts"]["p2830_full_carrier_separation_witness"])
        self.assertTrue(matrix["facts"]["graph_index_coverage_exact"])
        self.assertTrue(matrix["facts"]["obligation_matrix_executed"])
        self.assertFalse(matrix["facts"]["some_candidate_accepted"])
        self.assertFalse(matrix["facts"]["strict_graph_source_law_exported"])
        self.assertFalse(matrix["facts"]["typed_graph_to_K_or_Ltotal_coupling_exported"])
        self.assertFalse(matrix["facts"]["units_and_variational_derivative_exported"])
        self.assertFalse(matrix["accepted_as_source_law_coupling_theorem"])
        self.assertTrue(matrix["accepted_as_bounded_obligation_no_go"])

    def test_documents_updated(self):
        self.assertIn("P2831/S1781", MD_PATH.read_text(encoding="utf-8"))
        self.assertIn("P2831/S1781", (ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md").read_text(encoding="utf-8"))
        self.assertIn("P2831/S1781", (ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md").read_text(encoding="utf-8"))
        self.assertIn("P2831", (REPO / "AGENTS.md").read_text(encoding="utf-8"))


if __name__ == "__main__":
    unittest.main()
