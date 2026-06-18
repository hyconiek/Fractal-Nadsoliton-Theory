import json
import subprocess
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
SCRIPT = ROOT / "p2834_s1784_residual_second_variation_two_edge_toggle_audit.py"
JSON_PATH = ROOT / "generated" / "p2834_s1784_residual_second_variation_two_edge_toggle_audit.json"
MD_PATH = ROOT / "generated" / "p2834_s1784_residual_second_variation_two_edge_toggle_audit.md"
MANIFEST_PATH = ROOT / "generated" / "p2834_s1784_residual_second_variation_digest_manifest.json"


class P2834ResidualSecondVariationTwoEdgeToggleAuditTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        subprocess.run([sys.executable, str(SCRIPT)], check=True, cwd=ROOT)
        cls.payload = json.loads(JSON_PATH.read_text(encoding="utf-8"))
        cls.manifest = json.loads(MANIFEST_PATH.read_text(encoding="utf-8"))

    def test_status_and_inputs(self):
        self.assertEqual(
            self.payload["status"],
            "P2834_RESIDUAL_SECOND_VARIATION_TWO_EDGE_TOGGLE_WITNESS_NO_COUPLING_NO_CLOSURE",
        )
        self.assertEqual(
            self.payload["input_statuses"]["P2832"],
            "P2832_COMMON_NEIGHBOR_CYCLE_SOURCE_FORMULA_NO_GO_NO_COUPLING_NO_CLOSURE",
        )
        self.assertEqual(
            self.payload["input_statuses"]["P2833"],
            "P2833_EDGE_TOGGLE_RESPONSE_POLYNOMIAL_RESIDUAL_NO_GO_NO_COUPLING_NO_CLOSURE",
        )

    def test_residual_second_variation_counts(self):
        audit = self.payload["residual_second_variation_two_edge_toggle_audit"]
        self.assertEqual(audit["decoded_full_carrier_graph_count"], 16828)
        self.assertEqual(audit["p2833_residual_class_count_recomputed"], 67)
        self.assertEqual(audit["p2833_residual_graph_count_recomputed"], 138)
        self.assertTrue(audit["residual_only_audit"])
        self.assertFalse(audit["full_carrier_replay_performed"])
        self.assertTrue(audit["non_label_second_variation_formula_exported_for_candidate"])
        self.assertEqual(audit["two_edge_toggle_pair_count_per_graph"], 7140)
        self.assertEqual(audit["refined_residual_class_count"], 138)
        self.assertEqual(audit["refined_residual_collision_class_count"], 0)
        self.assertEqual(audit["refined_residual_collision_graph_count"], 0)
        self.assertEqual(audit["refined_residual_max_class_size"], 1)
        self.assertEqual(audit["refined_residual_defect_after_formula"], 0)
        self.assertEqual(audit["refined_class_size_histogram"]["1"], 138)

    def test_manifest_shape(self):
        audit = self.payload["residual_second_variation_two_edge_toggle_audit"]
        self.assertEqual(audit["manifest_path"], "fundamental_action_reconstruction/generated/p2834_s1784_residual_second_variation_digest_manifest.json")
        self.assertEqual(self.manifest["status"], "P2834_RESIDUAL_SECOND_VARIATION_MANIFEST")
        self.assertEqual(self.manifest["p2833_residual_class_count"], 67)
        self.assertEqual(self.manifest["p2833_residual_graph_count"], 138)
        self.assertEqual(self.manifest["two_edge_toggle_pair_count_per_graph"], 7140)
        self.assertEqual(len(self.manifest["graph_profile_digests"]), 138)
        self.assertEqual(len({row["second_variation_profile_sha256"] for row in self.manifest["graph_profile_digests"]}), 138)
        self.assertTrue(all(row["profile_row_count"] > 0 for row in self.manifest["graph_profile_digests"]))

    def test_acceptance_boundaries(self):
        matrix = self.payload["acceptance_matrix"]
        facts = matrix["facts"]
        self.assertTrue(facts["p2833_residual_boundary_reused"])
        self.assertTrue(facts["residual_only_audit"])
        self.assertFalse(facts["full_carrier_replay_performed"])
        self.assertTrue(facts["non_label_second_variation_formula_exported_for_candidate"])
        self.assertTrue(facts["p2833_residual_counts_match"])
        self.assertTrue(facts["second_variation_separates_p2833_residuals"])
        self.assertFalse(facts["full_carrier_source_law_coupling_exported"])
        self.assertFalse(facts["proved_variational_derivative_into_K_or_Ltotal_exported"])
        self.assertFalse(facts["typed_graph_to_K_or_Ltotal_coupling_theorem_exported"])
        self.assertFalse(facts["selector_bridge_or_role_transfer_imported"])
        self.assertTrue(matrix["accepted_as_residual_second_variation_witness"])
        self.assertFalse(matrix["accepted_as_source_law_coupling"])
        self.assertTrue(matrix["accepted_as_no_coupling_boundary"])

    def test_documents_updated(self):
        self.assertIn("P2834/S1784", MD_PATH.read_text(encoding="utf-8"))
        self.assertIn("P2834/S1784", (ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md").read_text(encoding="utf-8"))
        self.assertIn("P2834/S1784", (ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md").read_text(encoding="utf-8"))
        self.assertIn("P2834", (REPO / "AGENTS.md").read_text(encoding="utf-8"))


if __name__ == "__main__":
    unittest.main()
