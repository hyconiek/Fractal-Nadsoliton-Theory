import json
import subprocess
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
P2872_SCRIPT = ROOT / "p2872_s1822_cyclic_predecessor_endpoint_orientation_no_source_audit.py"
SCRIPT = ROOT / "p2873_s1823_exhaustive_dihedral_z12_endpoint_predicate_no_go_audit.py"
JSON_PATH = ROOT / "generated" / "p2873_s1823_exhaustive_dihedral_z12_endpoint_predicate_no_go_audit.json"
MD_PATH = ROOT / "generated" / "p2873_s1823_exhaustive_dihedral_z12_endpoint_predicate_no_go_audit.md"


class P2873ExhaustiveDihedralZ12EndpointPredicateNoGoAuditTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        for script in (P2872_SCRIPT, SCRIPT):
            subprocess.run([sys.executable, str(script)], check=True, cwd=ROOT, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
        cls.payload = json.loads(JSON_PATH.read_text(encoding="utf-8"))
        cls.audit = cls.payload["exhaustive_dihedral_z12_endpoint_predicate_no_go_audit"]

    def test_status_and_p2872_input(self):
        self.assertEqual(
            self.payload["status"],
            "P2873_EXHAUSTIVE_DIHEDRAL_Z12_ENDPOINT_PREDICATE_NO_GO_AUDIT_NO_CLOSURE",
        )
        self.assertEqual(
            self.audit["input_status_rechecked"],
            "P2872_CYCLIC_PREDECESSOR_ENDPOINT_ORIENTATION_NO_SOURCE_AUDIT_NO_CLOSURE",
        )

    def test_all_endpoint_predicates_are_enumerated(self):
        self.assertEqual(self.audit["predicate_count"], 4096)
        self.assertTrue(self.payload["acceptance_matrix"]["facts"]["all_4096_boolean_endpoint_predicates_enumerated"])

    def test_full_dihedral_invariants_cannot_select_endpoint(self):
        self.assertEqual(self.audit["full_dihedral_invariant_count"], 2)
        self.assertEqual(self.audit["full_dihedral_invariant_predicates"], [[], list(range(12))])
        facts = self.payload["acceptance_matrix"]["facts"]
        self.assertTrue(facts["full_dihedral_invariant_predicates_are_only_empty_and_all"])
        self.assertTrue(facts["singleton_11_not_full_dihedral_invariant"])
        self.assertTrue(facts["translation_orbit_of_11_is_all_z12"])

    def test_reflection_only_boundary_still_does_not_choose_11(self):
        self.assertEqual(self.audit["reflection0_invariant_count"], 128)
        self.assertIn([1, 11], self.audit["reflection0_atoms"])
        self.assertEqual(self.audit["singleton_11_reflection0_orbit"], [[1], [11]])
        self.assertEqual(self.audit["adjacent_pair_reflection0_orbit"], [[1, 11]])
        facts = self.payload["acceptance_matrix"]["facts"]
        self.assertTrue(facts["singleton_11_not_reflection0_invariant"])
        self.assertTrue(facts["adjacent_pair_reflection0_invariant_but_not_singleton"])

    def test_no_false_exports_and_documents_updated(self):
        flags = self.payload["decision"]["negative_export_flags"]
        self.assertFalse(flags["boundary_source_law_exported"])
        self.assertFalse(flags["orientation_source_law_exported"])
        self.assertFalse(flags["selector_or_localizer_source_exported"])
        self.assertFalse(flags["unit_bearing_coupling_localization_theorem_exported"])
        self.assertFalse(flags["ltotal_exported"])
        self.assertFalse(flags["toe_closure_exported"])
        self.assertIn("P2873/S1823", MD_PATH.read_text(encoding="utf-8"))
        self.assertIn("P2873/S1823", (ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md").read_text(encoding="utf-8"))
        self.assertIn("P2873/S1823", (ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md").read_text(encoding="utf-8"))
        self.assertIn("P2873", (REPO / "AGENTS.md").read_text(encoding="utf-8"))


if __name__ == "__main__":
    unittest.main()
