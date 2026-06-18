import json
import subprocess
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
P2870_SCRIPT = ROOT / "p2870_s1820_intrinsic_character_projector_selector_no_go_audit.py"
SCRIPT = ROOT / "p2871_s1821_exhaustive_gl22_invariant_projector_predicate_no_go_audit.py"
JSON_PATH = ROOT / "generated" / "p2871_s1821_exhaustive_gl22_invariant_projector_predicate_no_go_audit.json"
MD_PATH = ROOT / "generated" / "p2871_s1821_exhaustive_gl22_invariant_projector_predicate_no_go_audit.md"


class P2871ExhaustiveGL22InvariantProjectorPredicateNoGoAuditTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        for script in (P2870_SCRIPT, SCRIPT):
            subprocess.run([sys.executable, str(script)], check=True, cwd=ROOT, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
        cls.payload = json.loads(JSON_PATH.read_text(encoding="utf-8"))
        cls.audit = cls.payload["exhaustive_gl22_invariant_projector_predicate_no_go_audit"]

    def test_status_and_p2870_input(self):
        self.assertEqual(
            self.payload["status"],
            "P2871_EXHAUSTIVE_GL22_INVARIANT_PROJECTOR_PREDICATE_NO_GO_AUDIT_NO_CLOSURE",
        )
        self.assertEqual(
            self.audit["input_status_rechecked"],
            "P2870_INTRINSIC_CHARACTER_PROJECTOR_SELECTOR_NO_GO_AUDIT_NO_CLOSURE",
        )

    def test_exhausts_all_boolean_predicates(self):
        self.assertEqual(self.audit["predicate_count"], 16)
        facts = self.payload["acceptance_matrix"]["facts"]
        self.assertTrue(facts["all_16_boolean_predicates_enumerated"])

    def test_invariant_predicates_are_only_unions_of_identity_and_nonidentity_orbit(self):
        invariant_selected = sorted(record["selected"] for record in self.audit["invariant_predicates"])
        self.assertEqual(invariant_selected, [[], [1], [1, 5, 7, 11], [5, 7, 11]])
        self.assertEqual(len(self.audit["invariant_predicates"]), 4)

    def test_singleton_11_is_noninvariant_endpoint_import(self):
        target = self.audit["target_singleton_record"]
        self.assertEqual(target["selected"], [11])
        self.assertFalse(target["is_gl22_invariant"])
        self.assertEqual(target["orbit"], [[5], [7], [11]])
        facts = self.payload["acceptance_matrix"]["facts"]
        self.assertTrue(facts["no_invariant_predicate_selects_singleton_11"])
        self.assertEqual(self.audit["accepted_candidate_count"], 0)

    def test_no_false_exports_and_documents_updated(self):
        flags = self.payload["decision"]["negative_export_flags"]
        self.assertFalse(flags["boundary_source_law_exported"])
        self.assertFalse(flags["selector_or_localizer_source_exported"])
        self.assertFalse(flags["ltotal_exported"])
        self.assertFalse(flags["eom_closure_exported"])
        self.assertFalse(flags["toe_closure_exported"])
        self.assertIn("P2871/S1821", MD_PATH.read_text(encoding="utf-8"))
        self.assertIn("P2871/S1821", (ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md").read_text(encoding="utf-8"))
        self.assertIn("P2871/S1821", (ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md").read_text(encoding="utf-8"))
        self.assertIn("P2871", (REPO / "AGENTS.md").read_text(encoding="utf-8"))


if __name__ == "__main__":
    unittest.main()
