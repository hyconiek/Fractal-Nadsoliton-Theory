import json
import subprocess
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
SCRIPT = ROOT / "p2749_s1699_minimal_inversion_odd_source_coupling_polarity_audit.py"
OUT = ROOT / "generated" / "p2749_s1699_minimal_inversion_odd_source_coupling_polarity_audit.json"
MD = ROOT / "generated" / "p2749_s1699_minimal_inversion_odd_source_coupling_polarity_audit.md"


class P2749MinimalInversionOddSourceCouplingPolarityAuditTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        subprocess.run([sys.executable, str(SCRIPT)], cwd=ROOT, check=True)
        cls.payload = json.loads(OUT.read_text(encoding="utf-8"))
        cls.scan = cls.payload["content_evidence_scan"]
        cls.audit = cls.payload["coupling_audit"]
        cls.acceptance = cls.payload["acceptance_matrix"]

    def test_content_first_grep_detects_minimal_odd_source_boundaries(self):
        self.assertEqual(self.payload["status"], "P2749_MINIMAL_ODD_SOURCE_COUPLING_POLARITY_GAP")
        self.assertEqual(self.scan["content_pattern_count"], 4)
        self.assertTrue(self.scan["all_patterns_have_hits"])
        self.assertGreater(self.scan["hit_counts"]["post_p2748_missing_source"], 0)
        self.assertGreater(self.scan["hit_counts"]["coupling_polarity_boundary"], 0)

    def test_minimal_odd_source_has_two_equivariant_maps(self):
        self.assertEqual(self.audit["orientation_reversing_units"], [7, 11])
        self.assertEqual(self.audit["all_set_maps_count"], 4)
        self.assertEqual(self.audit["equivariant_map_count"], 2)
        self.assertEqual(self.audit["equivariant_maps"], [{"-1": -1, "1": 1}, {"-1": 1, "1": -1}])

    def test_p2721_polarity_exchanges_couplings(self):
        self.assertEqual(len(self.audit["p2721_polarity_rows"]), 4)
        self.assertEqual(self.audit["unique_coupled_map_count_after_p2721_polarity"], 2)
        self.assertIn("exchanges the two equivariant maps", self.audit["polarity_pairing_witness"])

    def test_acceptance_blocks_abstract_odd_source_export(self):
        self.assertFalse(self.acceptance["accepted_as_selector_source"])
        self.assertTrue(self.acceptance["facts"]["minimal_odd_source_representation_admitted"])
        self.assertFalse(self.acceptance["facts"]["coupling_polarity_unique"])
        self.assertIn("strict_source_sign_value_exported", self.acceptance["missing_criteria"])
        self.assertIn("p2721_polarity_coupling_exported", self.acceptance["missing_criteria"])
        self.assertTrue(all(value is False for value in self.payload["decision"]["negative_export_flags"].values()))

    def test_recommendation_and_docs_are_updated(self):
        self.assertIn("concrete strict source sign value", self.payload["decision"]["next_honest_step"])
        self.assertIn("P2697-P2749", self.payload["decision"]["next_honest_step"])
        self.assertIn("P2749/S1699", MD.read_text(encoding="utf-8"))
        self.assertIn("P2749/S1699", (ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md").read_text(encoding="utf-8"))
        self.assertIn("P2749/S1699", (ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md").read_text(encoding="utf-8"))
        self.assertIn("P2749", (REPO / "AGENTS.md").read_text(encoding="utf-8"))


if __name__ == "__main__":
    unittest.main()
