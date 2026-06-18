import json
import subprocess
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
P2869_SCRIPT = ROOT / "p2869_s1819_aut_character_idempotent_endpoint_localizer_no_source_audit.py"
SCRIPT = ROOT / "p2870_s1820_intrinsic_character_projector_selector_no_go_audit.py"
JSON_PATH = ROOT / "generated" / "p2870_s1820_intrinsic_character_projector_selector_no_go_audit.json"
MD_PATH = ROOT / "generated" / "p2870_s1820_intrinsic_character_projector_selector_no_go_audit.md"


class P2870IntrinsicCharacterProjectorSelectorNoGoAuditTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        for script in (P2869_SCRIPT, SCRIPT):
            subprocess.run([sys.executable, str(script)], check=True, cwd=ROOT, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
        cls.payload = json.loads(JSON_PATH.read_text(encoding="utf-8"))
        cls.audit = cls.payload["intrinsic_character_projector_selector_no_go_audit"]

    def test_status_and_p2869_input(self):
        self.assertEqual(
            self.payload["status"],
            "P2870_INTRINSIC_CHARACTER_PROJECTOR_SELECTOR_NO_GO_AUDIT_NO_CLOSURE",
        )
        self.assertEqual(
            self.audit["input_status_rechecked"],
            "P2869_AUT_CHARACTER_IDEMPOTENT_ENDPOINT_LOCALIZER_NO_SOURCE_AUDIT_NO_CLOSURE",
        )

    def test_v4_relabeling_ties_nonidentity_projectors(self):
        self.assertEqual(self.audit["nonidentity_projector_orbit"], [5, 7, 11])
        self.assertEqual(len(self.audit["v4_relabelings"]), 6)

    def test_intrinsic_features_do_not_single_out_11(self):
        projectors = self.audit["projector_table"]
        f5 = projectors["5"]["intrinsic_features"]
        f7 = projectors["7"]["intrinsic_features"]
        f11 = projectors["11"]["intrinsic_features"]
        for key in ("support_size", "positive_count", "negative_count", "zero_count", "coefficient_multiset", "l2_norm_squared", "scaled_denominator_support"):
            self.assertEqual(f5[key], f7[key])
            self.assertEqual(f7[key], f11[key])
        self.assertEqual(f5["moments"], f7["moments"])
        self.assertEqual(f7["moments"], f11["moments"])

    def test_target_match_is_only_named_boundary_target(self):
        self.assertEqual(self.audit["target_prime_vector_hits"], [11])
        facts = self.payload["acceptance_matrix"]["facts"]
        self.assertTrue(facts["target_match_singles_11_only_by_reusing_named_target"])
        self.assertTrue(facts["intrinsic_features_do_not_single_out_11_among_nonidentity_projectors"])
        self.assertEqual(self.audit["accepted_candidate_count"], 0)

    def test_no_false_exports_and_documents_updated(self):
        flags = self.payload["decision"]["negative_export_flags"]
        self.assertFalse(flags["boundary_source_law_exported"])
        self.assertFalse(flags["selector_or_localizer_source_exported"])
        self.assertFalse(flags["ltotal_exported"])
        self.assertFalse(flags["eom_closure_exported"])
        self.assertFalse(flags["toe_closure_exported"])
        self.assertIn("P2870/S1820", MD_PATH.read_text(encoding="utf-8"))
        self.assertIn("P2870/S1820", (ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md").read_text(encoding="utf-8"))
        self.assertIn("P2870/S1820", (ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md").read_text(encoding="utf-8"))
        self.assertIn("P2870", (REPO / "AGENTS.md").read_text(encoding="utf-8"))


if __name__ == "__main__":
    unittest.main()
