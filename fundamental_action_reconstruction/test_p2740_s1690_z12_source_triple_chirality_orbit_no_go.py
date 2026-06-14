import json
import subprocess
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
SCRIPT = ROOT / "p2740_s1690_z12_source_triple_chirality_orbit_no_go.py"
OUT = ROOT / "generated" / "p2740_s1690_z12_source_triple_chirality_orbit_no_go.json"
MD = ROOT / "generated" / "p2740_s1690_z12_source_triple_chirality_orbit_no_go.md"


class P2740Z12SourceTripleChiralityOrbitNoGoTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        subprocess.run([sys.executable, str(SCRIPT)], cwd=ROOT, check=True)
        cls.payload = json.loads(OUT.read_text(encoding="utf-8"))
        cls.scan = cls.payload["content_evidence_scan"]
        cls.audit = cls.payload["triple_chirality_audit"]
        cls.acceptance = cls.payload["acceptance_matrix"]

    def test_content_first_grep_detects_boundaries(self):
        self.assertEqual(self.payload["status"], "P2740_Z12_SOURCE_TRIPLE_CHIRALITY_ORBIT_NO_GO")
        self.assertEqual(self.scan["content_pattern_count"], 4)
        self.assertTrue(self.scan["all_patterns_have_hits"])
        self.assertGreater(self.scan["hit_counts"]["existing_three_cycle_reference_boundary"], 0)
        self.assertGreater(self.scan["hit_counts"]["p2721_polarity_boundary"], 0)

    def test_ordered_triple_chirality_counts_are_balanced(self):
        self.assertEqual(self.audit["ordered_distinct_triples"], 1320)
        self.assertEqual(self.audit["unordered_triples"], 220)
        self.assertEqual(self.audit["positive_ordered_triples"], 660)
        self.assertEqual(self.audit["negative_ordered_triples"], 660)
        self.assertEqual(self.audit["translation_orbits_with_nonzero_signed_sum"], 0)
        self.assertEqual(self.audit["affine_orbits_with_nonzero_signed_sum"], 0)

    def test_orbit_structure_is_finite_and_nonempty(self):
        self.assertGreater(self.audit["translation_unordered_orbit_count"], 0)
        self.assertGreater(self.audit["affine_unordered_orbit_count"], 0)
        self.assertTrue(all(row["signed_sum"] == 0 for row in self.audit["translation_orbit_rows"]))
        self.assertTrue(all(row["signed_sum"] == 0 for row in self.audit["affine_orbit_rows"]))

    def test_acceptance_blocks_closure_and_source_export(self):
        self.assertFalse(self.acceptance["accepted_as_strict_signed_source"])
        self.assertIn("translation_orbit_signed_source_survives", self.acceptance["missing_criteria"])
        self.assertIn("strict_source_localizer_exported", self.acceptance["missing_criteria"])
        self.assertIn("p2721_polarity_coupling_exported", self.acceptance["missing_criteria"])
        self.assertTrue(all(value is False for value in self.payload["decision"]["negative_export_flags"].values()))

    def test_recommendation_and_docs_are_updated(self):
        self.assertIn("real nonzero pointwise sign", self.payload["decision"]["next_honest_step"])
        self.assertIn("P2697-P2740", self.payload["decision"]["next_honest_step"])
        self.assertIn("P2740/S1690", MD.read_text(encoding="utf-8"))
        self.assertIn("P2740/S1690", (ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md").read_text(encoding="utf-8"))
        self.assertIn("P2740/S1690", (ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md").read_text(encoding="utf-8"))
        self.assertIn("P2740", (REPO / "AGENTS.md").read_text(encoding="utf-8"))


if __name__ == "__main__":
    unittest.main()
