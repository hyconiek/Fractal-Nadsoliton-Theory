import json
import math
import subprocess
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
SCRIPT = ROOT / "p2839_s1789_graph_bit_localization_map_obstruction_audit.py"
JSON_PATH = ROOT / "generated" / "p2839_s1789_graph_bit_localization_map_obstruction_audit.json"
MD_PATH = ROOT / "generated" / "p2839_s1789_graph_bit_localization_map_obstruction_audit.md"


class P2839GraphBitLocalizationMapObstructionAuditTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        subprocess.run([sys.executable, str(SCRIPT)], check=True, cwd=ROOT)
        cls.payload = json.loads(JSON_PATH.read_text(encoding="utf-8"))

    def test_status(self):
        self.assertEqual(
            self.payload["status"],
            "P2839_LOCALIZATION_PULLBACK_OBSTRUCTION_NO_GO_NO_COUPLING_NO_CLOSURE",
        )

    def test_label_gauge_inventory(self):
        audit = self.payload["graph_bit_localization_map_obstruction_audit"]
        inventory = audit["finite_label_gauge_inventory"]
        self.assertEqual(inventory["vertex_count"], 16)
        self.assertEqual(inventory["edge_bit_count"], 120)
        self.assertEqual(inventory["label_permutation_count_exact"], math.factorial(16))
        self.assertTrue(inventory["decoded_graphs_are_finite_unlabeled_isomorphism_class_representatives"])
        self.assertFalse(inventory["canonical_vertex_coordinate_source_exported"])

    def test_candidate_localization_maps(self):
        audit = self.payload["graph_bit_localization_map_obstruction_audit"]
        candidates = {row["candidate"]: row for row in audit["candidate_localization_maps"]}
        self.assertEqual(set(candidates), {
            "label_index_to_fixed_points",
            "edge_orbit_to_local_source_bins",
            "graphon_step_function",
            "spectral_embedding_localization",
        })
        self.assertTrue(all(row["premises"]["graph_bit_domain"] for row in candidates.values()))
        self.assertFalse(any(row["accepted_as_localization_map"] for row in candidates.values()))
        self.assertFalse(candidates["label_index_to_fixed_points"]["premises"]["label_gauge_invariant_vertex_or_edge_coordinates"])
        self.assertTrue(candidates["edge_orbit_to_local_source_bins"]["premises"]["label_gauge_invariant_vertex_or_edge_coordinates"])
        self.assertFalse(candidates["edge_orbit_to_local_source_bins"]["premises"]["target_field_or_spacetime_support"])

    def test_obstruction_summary(self):
        audit = self.payload["graph_bit_localization_map_obstruction_audit"]
        summary = audit["localization_obstruction_summary"]
        self.assertEqual(summary["candidate_count"], 4)
        self.assertEqual(summary["accepted_localization_map_count"], 0)
        self.assertEqual(summary["accepted_localization_maps"], [])
        self.assertIn("locality_covariance_rule", summary["common_hard_blockers"])
        self.assertIn("target_independent_units_or_measure", summary["common_hard_blockers"])
        self.assertIn("compatibility_with_variational_chain_rule", summary["common_hard_blockers"])
        self.assertEqual(summary["blocker_histogram"]["label_gauge_invariant_vertex_or_edge_coordinates"], 3)
        self.assertEqual(summary["blocker_histogram"]["pullback_formula_Aij_to_field_data"], 2)

    def test_acceptance_boundaries(self):
        matrix = self.payload["acceptance_matrix"]
        facts = matrix["facts"]
        self.assertTrue(facts["p2835_combined_separator_rechecked"])
        self.assertTrue(facts["at_least_one_localization_candidate_tested"])
        self.assertFalse(facts["accepted_localization_map_exported"])
        self.assertFalse(facts["canonical_vertex_coordinate_source_exported"])
        self.assertFalse(facts["selector_bridge_or_role_transfer_imported"])
        self.assertTrue(matrix["accepted_as_localization_obstruction_audit"])
        self.assertFalse(matrix["accepted_as_evaluation_pullback_localization_map"])
        self.assertTrue(matrix["accepted_as_localization_no_go"])

    def test_documents_updated(self):
        self.assertIn("P2839/S1789", MD_PATH.read_text(encoding="utf-8"))
        self.assertIn("P2839/S1789", (ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md").read_text(encoding="utf-8"))
        self.assertIn("P2839/S1789", (ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md").read_text(encoding="utf-8"))
        self.assertIn("P2839", (REPO / "AGENTS.md").read_text(encoding="utf-8"))


if __name__ == "__main__":
    unittest.main()
