import json
import math
import subprocess
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
SCRIPT = ROOT / "p2754_s1704_shannon_entropy_four_bit_selector_audit.py"
OUT = ROOT / "generated" / "p2754_s1704_shannon_entropy_four_bit_selector_audit.json"
MD = ROOT / "generated" / "p2754_s1704_shannon_entropy_four_bit_selector_audit.md"


class P2754ShannonEntropyFourBitSelectorAuditTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        subprocess.run([sys.executable, str(SCRIPT)], cwd=ROOT, check=True)
        cls.payload = json.loads(OUT.read_text(encoding="utf-8"))
        cls.scan = cls.payload["content_evidence_scan"]
        cls.audit = cls.payload["shannon_entropy_selector_audit"]
        cls.acceptance = cls.payload["acceptance_matrix"]

    def test_content_first_grep_detects_entropy_and_selector_boundaries(self):
        self.assertEqual(self.payload["status"], "P2754_SHANNON_ENTROPY_FOUR_BIT_SELECTOR_AUDIT_NO_GO")
        self.assertEqual(self.scan["content_pattern_count"], 4)
        self.assertTrue(self.scan["all_patterns_have_hits"])
        self.assertGreater(self.scan["hit_counts"]["post_p2753_pivot_boundary"], 0)
        self.assertGreater(self.scan["hit_counts"]["entropy_alpha_geo_lane"], 0)
        self.assertGreater(self.scan["hit_counts"]["selector_boundary"], 0)

    def test_four_ln2_is_verified_as_four_bit_entropy(self):
        self.assertTrue(self.audit["four_bit_entropy_matches_4_ln2"])
        self.assertAlmostEqual(self.audit["four_bit_entropy_nats"], 4 * math.log(2), places=12)
        self.assertAlmostEqual(self.audit["uniform_16_entropy_nats"], 4 * math.log(2), places=12)

    def test_z12_entropy_scan_is_inversion_invariant(self):
        scan = self.audit["z12_integer_weight_entropy_scan"]
        self.assertEqual(scan["quanta"], 4)
        self.assertEqual(scan["composition_count"], 1365)
        self.assertEqual(scan["inversion_entropy_failure_count"], 0)
        self.assertGreater(scan["distinct_entropy_value_count"], 1)
        self.assertAlmostEqual(scan["max_entropy_for_four_quanta_nats"], math.log(4), places=12)

    def test_entropy_has_zero_equivariant_maps_to_orientation_torsor(self):
        torsor = self.audit["torsor_equivariance_test"]
        self.assertEqual(torsor["orientation_reversing_units"], [7, 11])
        self.assertEqual(torsor["torsor_fixed_points_under_reversing_units"], 0)
        self.assertEqual(torsor["equivariant_maps_from_four_bit_max_entropy_singleton"], 0)
        self.assertEqual(torsor["equivariant_maps_from_entropy_value_quotient"], 0)

    def test_acceptance_blocks_scalar_entropy_selector_promotion(self):
        self.assertFalse(self.acceptance["accepted_as_entropy_generated_selector"])
        self.assertTrue(self.acceptance["facts"]["four_ln2_is_verified_as_four_bit_shannon_entropy"])
        self.assertTrue(self.acceptance["facts"]["z12_entropy_is_inversion_invariant_on_integer_weight_scan"])
        self.assertEqual(self.acceptance["facts"]["equivariant_maps_from_entropy_to_orientation_torsor"], 0)
        self.assertFalse(self.acceptance["facts"]["new_inversion_odd_entropy_current_exported"])
        self.assertFalse(self.acceptance["facts"]["p2721_entropy_coupling_theorem_exported"])
        self.assertTrue(all(value is False for value in self.payload["decision"]["negative_export_flags"].values()))

    def test_recommendation_and_docs_are_updated(self):
        self.assertIn("P2697-P2754", self.payload["decision"]["next_honest_step"])
        self.assertIn("P2754/S1704", MD.read_text(encoding="utf-8"))
        self.assertIn("P2754/S1704", (ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md").read_text(encoding="utf-8"))
        self.assertIn("P2754/S1704", (ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md").read_text(encoding="utf-8"))
        self.assertIn("P2754", (REPO / "AGENTS.md").read_text(encoding="utf-8"))


if __name__ == "__main__":
    unittest.main()
