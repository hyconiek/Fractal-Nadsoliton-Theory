import json
import subprocess
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
SCRIPT = ROOT / "p3146_s2096_axiom_unit_action_measure_bridge.py"
OUT = ROOT / "generated" / "p3146_s2096_axiom_unit_action_measure_bridge.json"
MD = ROOT / "generated" / "p3146_s2096_axiom_unit_action_measure_bridge.md"


class P3146AxiomUnitActionMeasureBridgeTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        subprocess.run([sys.executable, str(SCRIPT)], cwd=ROOT, check=True, stdout=subprocess.PIPE, text=True)
        cls.payload = json.loads(OUT.read_text(encoding="utf-8"))

    def test_status_and_counts(self):
        self.assertEqual(self.payload["status"], "P3146_AXIOM_UNIT_ACTION_MEASURE_BRIDGE_CONDITIONAL_NON_STRICT")
        self.assertTrue(self.payload["input_hashes"]["P3145"])
        counts = self.payload["finite_theorem"]["finite_counts"]
        self.assertEqual(counts["axiom_subsets_audited"], 8)
        self.assertEqual(counts["dimension_basis_size"], 3)
        self.assertEqual(counts["fully_covering_subsets"], 1)
        self.assertEqual(counts["minimal_fully_covering_subsets"], 1)
        self.assertEqual(counts["strict_source_subsets"], 0)
        self.assertEqual(counts["conditional_unit_action_measure_subsets"], 1)

    def test_only_full_triple_covers_unit_targets(self):
        full = [row for row in self.payload["audit_rows"] if row["all_unit_targets_covered"]]
        self.assertEqual(len(full), 1)
        self.assertEqual(set(full[0]["axiom_set"]), {"A_cell", "A_clock", "A_action"})
        self.assertEqual(full[0]["rank_HTL"], 3)
        self.assertTrue(full[0]["unit_bearing_action_measure_available_conditionally"])
        self.assertFalse(full[0]["strict_source_exported"])

    def test_decision_preserves_boundaries(self):
        decision = self.payload["decision"]
        self.assertIn("full triple", decision["bounded_result"])
        self.assertIn("imported postulates", decision["why_not_strict"])
        self.assertIn("A_action", decision["next_honest_step"])
        self.assertTrue(all(value is False for value in decision["negative_export_flags"].values()))

    def test_docs_updated(self):
        self.assertIn("P3146/S2096", MD.read_text(encoding="utf-8"))
        self.assertIn("P3146/S2096", (ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md").read_text(encoding="utf-8"))
        self.assertIn("P3146/S2096", (ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md").read_text(encoding="utf-8"))
        self.assertIn("P3146", (REPO / "AGENTS.md").read_text(encoding="utf-8"))


if __name__ == "__main__":
    unittest.main()
