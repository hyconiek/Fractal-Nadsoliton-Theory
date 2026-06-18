import json
import subprocess
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
P2871_SCRIPT = ROOT / "p2871_s1821_exhaustive_gl22_invariant_projector_predicate_no_go_audit.py"
SCRIPT = ROOT / "p2872_s1822_cyclic_predecessor_endpoint_orientation_no_source_audit.py"
JSON_PATH = ROOT / "generated" / "p2872_s1822_cyclic_predecessor_endpoint_orientation_no_source_audit.json"
MD_PATH = ROOT / "generated" / "p2872_s1822_cyclic_predecessor_endpoint_orientation_no_source_audit.md"


class P2872CyclicPredecessorEndpointOrientationNoSourceAuditTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        for script in (P2871_SCRIPT, SCRIPT):
            subprocess.run([sys.executable, str(script)], check=True, cwd=ROOT, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
        cls.payload = json.loads(JSON_PATH.read_text(encoding="utf-8"))
        cls.audit = cls.payload["cyclic_predecessor_endpoint_orientation_no_source_audit"]

    def test_status_and_p2871_input(self):
        self.assertEqual(
            self.payload["status"],
            "P2872_CYCLIC_PREDECESSOR_ENDPOINT_ORIENTATION_NO_SOURCE_AUDIT_NO_CLOSURE",
        )
        self.assertEqual(
            self.audit["input_status_rechecked"],
            "P2871_EXHAUSTIVE_GL22_INVARIANT_PROJECTOR_PREDICATE_NO_GO_AUDIT_NO_CLOSURE",
        )

    def test_predecessor_11_exactly_represents_target(self):
        self.assertEqual(self.audit["predecessor_endpoint"], 11)
        self.assertTrue(self.audit["exact_predecessor_target_representation"])
        vector = self.audit["predecessor_prime_vector_scaled_by_9_over_5"]
        self.assertEqual(vector["11"]["fraction"], "9/5")
        self.assertEqual(vector["5"]["fraction"], "0/1")
        self.assertEqual(vector["7"]["fraction"], "0/1")

    def test_reflection_invariant_adjacent_predicates_cannot_select_11(self):
        invariant_selected = sorted(record["selected"] for record in self.audit["reflection_invariant_adjacent_predicates"])
        self.assertEqual(invariant_selected, [[], [1, 11]])
        facts = self.payload["acceptance_matrix"]["facts"]
        self.assertTrue(facts["reflection_swaps_successor_1_and_predecessor_11"])
        self.assertTrue(facts["reflection_invariant_predicates_do_not_select_singleton_11"])
        self.assertEqual(self.audit["accepted_candidate_count"], 0)

    def test_successor_endpoint_is_target_blind(self):
        vector = self.audit["successor_prime_vector_scaled_by_9_over_5"]
        self.assertTrue(all(value["fraction"] == "0/1" for value in vector.values()))
        self.assertTrue(self.payload["acceptance_matrix"]["facts"]["successor_endpoint_is_target_blind"])

    def test_no_false_exports_and_documents_updated(self):
        flags = self.payload["decision"]["negative_export_flags"]
        self.assertFalse(flags["boundary_source_law_exported"])
        self.assertFalse(flags["orientation_source_law_exported"])
        self.assertFalse(flags["selector_or_localizer_source_exported"])
        self.assertFalse(flags["ltotal_exported"])
        self.assertFalse(flags["toe_closure_exported"])
        self.assertIn("P2872/S1822", MD_PATH.read_text(encoding="utf-8"))
        self.assertIn("P2872/S1822", (ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md").read_text(encoding="utf-8"))
        self.assertIn("P2872/S1822", (ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md").read_text(encoding="utf-8"))
        self.assertIn("P2872", (REPO / "AGENTS.md").read_text(encoding="utf-8"))


if __name__ == "__main__":
    unittest.main()
