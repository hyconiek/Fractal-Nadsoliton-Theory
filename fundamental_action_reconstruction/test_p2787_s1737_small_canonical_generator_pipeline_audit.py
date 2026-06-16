import json
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
GEN = ROOT / "generated"
JSON = GEN / "p2787_s1737_small_canonical_generator_pipeline_audit.json"
MD = GEN / "p2787_s1737_small_canonical_generator_pipeline_audit.md"


class P2787SmallCanonicalGeneratorPipelineAuditTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        if not JSON.exists():
            import p2787_s1737_small_canonical_generator_pipeline_audit as p2787
            p2787.main()
        cls.payload = json.loads(JSON.read_text(encoding="utf-8"))

    def test_status_and_inputs(self):
        self.assertEqual(self.payload["status"], "P2787_SMALL_CANONICAL_GENERATOR_PIPELINE_AUDIT_NO_CLOSURE")
        self.assertEqual(self.payload["input_statuses"]["P2786"], "P2786_GRAPH6_PROVENANCE_TOOLCHAIN_GATE_NO_CLOSURE")
        self.assertIn("complete smaller regular graph class", self.payload["audited_question"])

    def test_complete_small_pipeline_counts(self):
        witness = self.payload["small_pipeline_witness"]
        self.assertEqual(witness["labeled_candidate_count"], 19355)
        self.assertEqual(witness["connected_labeled_candidate_count"], 19355)
        self.assertEqual(witness["isomorphism_class_count"], 6)
        self.assertEqual(witness["pair_count_after_quotient"], 15)
        self.assertEqual(witness["exact_charpoly_collision_counts"], {
            "adjacency_charpoly_coefficients": 0,
            "laplacian_charpoly_coefficients": 0,
            "signless_laplacian_charpoly_coefficients": 0,
        })
        self.assertTrue(witness["all_pairs_separated_by_all_three_exact_charpolys"])
        self.assertEqual(len(witness["representative_rows"]), 6)
        self.assertEqual(len(witness["pair_rows"]), 15)
        self.assertEqual(sum(row["labeled_member_count"] for row in witness["representative_rows"]), 19355)
        for row in witness["representative_rows"]:
            self.assertEqual(row["edge_count"], 16)
            self.assertTrue(row["graph6"].startswith("G"))
            self.assertEqual(len(row["adjacency_charpoly_coefficients"]), 9)

    def test_acceptance_blocks_16node_and_closure(self):
        acceptance = self.payload["acceptance_matrix"]
        self.assertTrue(acceptance["accepted_as_small_self_contained_generator_pipeline_certificate"])
        self.assertFalse(acceptance["accepted_as_full_16node_canonical_generator_certificate"])
        self.assertFalse(acceptance["accepted_as_strict_spectral_source_law"])
        self.assertFalse(acceptance["accepted_as_canonical_nadsoliton_geometry_source"])
        self.assertFalse(acceptance["accepted_as_ltotal_or_toe_promotion"])
        self.assertIn("P2697-P2787", self.payload["decision"]["next_honest_step"])

    def test_documentation_updated(self):
        self.assertIn("P2787/S1737", MD.read_text(encoding="utf-8"))
        self.assertIn("P2787/S1737", (ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md").read_text(encoding="utf-8"))
        self.assertIn("P2787/S1737", (ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md").read_text(encoding="utf-8"))
        self.assertIn("P2787", (REPO / "AGENTS.md").read_text(encoding="utf-8"))


if __name__ == "__main__":
    unittest.main()
