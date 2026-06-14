import json
import subprocess
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
SCRIPT = ROOT / "p2742_s1692_source_triple_affine_weighted_signed_aggregate_no_go.py"
OUT = ROOT / "generated" / "p2742_s1692_source_triple_affine_weighted_signed_aggregate_no_go.json"
MD = ROOT / "generated" / "p2742_s1692_source_triple_affine_weighted_signed_aggregate_no_go.md"


class P2742SourceTripleAffineWeightedSignedAggregateNoGoTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        subprocess.run([sys.executable, str(SCRIPT)], cwd=ROOT, check=True)
        cls.payload = json.loads(OUT.read_text(encoding="utf-8"))
        cls.scan = cls.payload["content_evidence_scan"]
        cls.audit = cls.payload["aggregate_audit"]
        cls.acceptance = cls.payload["acceptance_matrix"]

    def test_content_first_grep_detects_weighting_escape_hatch(self):
        self.assertEqual(self.payload["status"], "P2742_SOURCE_TRIPLE_AFFINE_WEIGHTED_SIGNED_AGGREGATE_NO_GO")
        self.assertEqual(self.scan["content_pattern_count"], 4)
        self.assertTrue(self.scan["all_patterns_have_hits"])
        self.assertGreater(self.scan["hit_counts"]["p2741_localizer_no_go"], 0)
        self.assertGreater(self.scan["hit_counts"]["affine_orbit_weighting_boundary"], 0)

    def test_affine_orbit_signed_coefficients_all_vanish(self):
        self.assertEqual(self.audit["ordered_distinct_triples"], 1320)
        self.assertEqual(self.audit["affine_group_size"], 48)
        self.assertEqual(self.audit["affine_ordered_orbit_count"], 34)
        self.assertEqual(self.audit["affine_orbit_sizes"], [24, 48])
        self.assertEqual(self.audit["orbits_with_nonzero_signed_sum_coefficient"], 0)
        self.assertTrue(all(row["signed_sum_coefficient"] == 0 for row in self.audit["orbit_rows"]))

    def test_linear_map_and_pairing_witness(self):
        self.assertEqual(self.audit["signed_sum_linear_map_rank"], 0)
        self.assertEqual(self.audit["signed_sum_linear_map_nullity_over_orbit_weights"], 34)
        self.assertEqual(self.audit["one_hot_weight_crosscheck_nonzero_count"], 0)
        self.assertTrue(self.audit["all_unit_11_pairing_witnesses_pass"])
        self.assertTrue(all(row["unit_11_pairs_chirality_opposites"] for row in self.audit["orbit_rows"]))

    def test_acceptance_blocks_orbit_safe_source_export(self):
        self.assertFalse(self.acceptance["accepted_as_orbit_safe_signed_source"])
        self.assertIn("affine_orbit_weighted_signed_aggregate_nonzero", self.acceptance["missing_criteria"])
        self.assertIn("strict_source_localizer_exported", self.acceptance["missing_criteria"])
        self.assertIn("p2721_polarity_coupling_exported", self.acceptance["missing_criteria"])
        self.assertTrue(all(value is False for value in self.payload["decision"]["negative_export_flags"].values()))

    def test_recommendation_and_docs_are_updated(self):
        self.assertIn("Do not continue the P2740/P2741 source-triple lane", self.payload["decision"]["next_honest_step"])
        self.assertIn("P2697-P2742", self.payload["decision"]["next_honest_step"])
        self.assertIn("P2742/S1692", MD.read_text(encoding="utf-8"))
        self.assertIn("P2742/S1692", (ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md").read_text(encoding="utf-8"))
        self.assertIn("P2742/S1692", (ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md").read_text(encoding="utf-8"))
        self.assertIn("P2742", (REPO / "AGENTS.md").read_text(encoding="utf-8"))


if __name__ == "__main__":
    unittest.main()
