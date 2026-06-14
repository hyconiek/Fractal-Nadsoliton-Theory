import json
import subprocess
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
SCRIPT = ROOT / "p2741_s1691_source_triple_localizer_fixed_point_no_go.py"
OUT = ROOT / "generated" / "p2741_s1691_source_triple_localizer_fixed_point_no_go.json"
MD = ROOT / "generated" / "p2741_s1691_source_triple_localizer_fixed_point_no_go.md"


class P2741SourceTripleLocalizerFixedPointNoGoTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        subprocess.run([sys.executable, str(SCRIPT)], cwd=ROOT, check=True)
        cls.payload = json.loads(OUT.read_text(encoding="utf-8"))
        cls.scan = cls.payload["content_evidence_scan"]
        cls.audit = cls.payload["localizer_audit"]
        cls.acceptance = cls.payload["acceptance_matrix"]

    def test_content_first_grep_detects_localizer_boundaries(self):
        self.assertEqual(self.payload["status"], "P2741_SOURCE_TRIPLE_LOCALIZER_FIXED_POINT_NO_GO")
        self.assertEqual(self.scan["content_pattern_count"], 4)
        self.assertTrue(self.scan["all_patterns_have_hits"])
        self.assertGreater(self.scan["hit_counts"]["p2740_localizer_obligation"], 0)
        self.assertGreater(self.scan["hit_counts"]["p2721_coupling_boundary"], 0)

    def test_fixed_point_orbit_counts_are_exact(self):
        self.assertEqual(self.audit["ordered_distinct_triples"], 1320)
        self.assertEqual(self.audit["translation_group_size"], 12)
        self.assertEqual(self.audit["affine_group_size"], 48)
        self.assertEqual(self.audit["translation_ordered_orbit_count"], 110)
        self.assertEqual(self.audit["translation_orbit_sizes"], [12])
        self.assertEqual(self.audit["translation_fixed_ordered_triples"], 0)
        self.assertEqual(self.audit["translation_singleton_orbits"], 0)

    def test_affine_has_no_singleton_localizer_orbit(self):
        self.assertEqual(self.audit["affine_ordered_orbit_count"], 34)
        self.assertEqual(self.audit["affine_fixed_ordered_triples"], 0)
        self.assertEqual(self.audit["affine_singleton_orbits"], 0)
        self.assertTrue(min(self.audit["affine_orbit_sizes"]) > 1)

    def test_acceptance_blocks_localizer_and_closure(self):
        self.assertFalse(self.acceptance["accepted_as_strict_source_localizer"])
        self.assertIn("translation_safe_fixed_ordered_triple_exists", self.acceptance["missing_criteria"])
        self.assertIn("affine_safe_singleton_ordered_orbit_exists", self.acceptance["missing_criteria"])
        self.assertIn("strict_source_localizer_exported", self.acceptance["missing_criteria"])
        self.assertTrue(all(value is False for value in self.payload["decision"]["negative_export_flags"].values()))

    def test_recommendation_and_docs_are_updated(self):
        self.assertIn("no translation/affine-safe fixed ordered triple", self.payload["decision"]["next_honest_step"])
        self.assertIn("P2697-P2741", self.payload["decision"]["next_honest_step"])
        self.assertIn("P2741/S1691", MD.read_text(encoding="utf-8"))
        self.assertIn("P2741/S1691", (ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md").read_text(encoding="utf-8"))
        self.assertIn("P2741/S1691", (ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md").read_text(encoding="utf-8"))
        self.assertIn("P2741", (REPO / "AGENTS.md").read_text(encoding="utf-8"))


if __name__ == "__main__":
    unittest.main()
