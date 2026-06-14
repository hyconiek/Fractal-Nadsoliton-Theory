from __future__ import annotations

import json
import subprocess
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
SCRIPT = ROOT / "p2708_s1658_z12_boundary_1_cocycle_selector_source_obstruction.py"
OUT = ROOT / "generated" / "p2708_s1658_z12_boundary_1_cocycle_selector_source_obstruction.json"
MD = ROOT / "generated" / "p2708_s1658_z12_boundary_1_cocycle_selector_source_obstruction.md"


class P2708Z12BoundaryCocycleSelectorObstructionTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        subprocess.run([sys.executable, str(SCRIPT)], cwd=ROOT, check=True)
        cls.payload = json.loads(OUT.read_text(encoding="utf-8"))

    def test_cohomology_computation(self):
        cohomology = self.payload["cohomology_certificate"]
        self.assertEqual(cohomology["rank_d0"], 11)
        self.assertEqual(cohomology["h1_dimension"], 1)
        self.assertEqual(cohomology["primitive_boundary_pairing"], 12)

    def test_aut_inversion_kills_canonical_nonzero_orientation(self):
        aut = self.payload["automorphism_action_certificate"]
        self.assertTrue(aut["rotation_fixed"])
        self.assertTrue(aut["inversion_sends_omega_to_minus_omega"])
        self.assertFalse(aut["aut_z12_invariant_nonzero_orientation_exists"])
        self.assertEqual(aut["invariant_subspace_dimension_under_inversion"], 0)

    def test_decision_preserves_no_strict_selector(self):
        self.assertEqual(self.payload["status"], "P2708_Z12_BOUNDARY_1_COCYCLE_OBSTRUCTION_NO_STRICT_SELECTOR")
        self.assertTrue(all(not row["strict_provider_obligation_met"] for row in self.payload["premise_vs_strict_table"]))
        decision = self.payload["decision"]
        self.assertFalse(decision["boundary_1_cocycle_exports_nonpremise_selector"])
        self.assertTrue(all(value is False for value in decision["negative_export_flags"].values()))
        self.assertIn("missing sign", decision["next_honest_step"])

    def test_outputs_and_docs_written(self):
        self.assertIn("P2708/S1658", MD.read_text(encoding="utf-8"))
        self.assertIn("P2708/S1658", (ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md").read_text(encoding="utf-8"))
        self.assertIn("P2708/S1658", (ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md").read_text(encoding="utf-8"))
        self.assertIn("P2708", (REPO / "AGENTS.md").read_text(encoding="utf-8"))


if __name__ == "__main__":
    unittest.main()
