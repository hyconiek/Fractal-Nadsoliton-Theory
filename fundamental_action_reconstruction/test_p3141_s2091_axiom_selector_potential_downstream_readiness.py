import json
import subprocess
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
SCRIPT = ROOT / "p3141_s2091_axiom_selector_potential_downstream_readiness.py"
OUT = ROOT / "generated" / "p3141_s2091_axiom_selector_potential_downstream_readiness.json"
MD = ROOT / "generated" / "p3141_s2091_axiom_selector_potential_downstream_readiness.md"


class P3141AxiomSelectorPotentialDownstreamReadinessTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        subprocess.run([sys.executable, str(SCRIPT)], cwd=ROOT, check=True, stdout=subprocess.PIPE, text=True)
        cls.payload = json.loads(OUT.read_text(encoding="utf-8"))

    def test_status_inputs_and_counts(self):
        self.assertEqual(self.payload["status"], "P3141_AXIOM_SELECTOR_POTENTIAL_DOWNSTREAM_READINESS_BOUNDED_NON_STRICT")
        self.assertTrue(all(self.payload["input_hashes"].values()))
        counts = self.payload["finite_theorem"]["finite_counts"]
        self.assertEqual(counts["pair_space_size"], 24)
        self.assertEqual(counts["weight_pairs_tested"], 9)
        self.assertEqual(counts["unique_minimizer_rows"], 9)
        self.assertEqual(counts["strict_weight_source_rows"], 0)
        self.assertEqual(counts["unit_bearing_scale_rows"], 0)

    def test_potential_rows(self):
        rows = self.payload["potential_rows"]
        self.assertEqual(len(rows), 9)
        self.assertTrue(all(row["unique_minimizer"] for row in rows))
        self.assertTrue(all(row["selected_pair"] == [0, 1] for row in rows))
        self.assertTrue(all(row["first_gap"] > 0 for row in rows))
        self.assertTrue(all(not row["strict_weight_source_exported"] for row in rows))

    def test_downstream_gates_and_exports(self):
        gates = {row["gate"]: row["passed"] for row in self.payload["downstream_gate_rows"]}
        self.assertTrue(gates["finite_unique_minimizer"])
        self.assertTrue(gates["non_strict_axiom_consistency"])
        self.assertFalse(gates["strict_weight_source"])
        self.assertFalse(gates["unit_bearing_scale"])
        self.assertFalse(gates["field_variational_lift"])
        self.assertFalse(gates["bridge_role_transfer_ToE"])
        self.assertTrue(all(value is False for value in self.payload["decision"]["negative_export_flags"].values()))

    def test_docs_updated(self):
        self.assertIn("P3141/S2091", MD.read_text(encoding="utf-8"))
        self.assertIn("P3141/S2091", (ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md").read_text(encoding="utf-8"))
        self.assertIn("P3141/S2091", (ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md").read_text(encoding="utf-8"))
        self.assertIn("P3141", (REPO / "AGENTS.md").read_text(encoding="utf-8"))


if __name__ == "__main__":
    unittest.main()
