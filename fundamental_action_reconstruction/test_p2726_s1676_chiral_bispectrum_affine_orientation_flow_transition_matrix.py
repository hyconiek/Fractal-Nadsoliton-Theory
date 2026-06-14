from __future__ import annotations

import json
import subprocess
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
SCRIPT = ROOT / "p2726_s1676_chiral_bispectrum_affine_orientation_flow_transition_matrix.py"
OUT = ROOT / "generated" / "p2726_s1676_chiral_bispectrum_affine_orientation_flow_transition_matrix.json"
MD = ROOT / "generated" / "p2726_s1676_chiral_bispectrum_affine_orientation_flow_transition_matrix.md"


class P2726ChiralBispectrumAffineOrientationFlowTransitionMatrixTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        subprocess.run([sys.executable, str(SCRIPT)], cwd=ROOT, check=True)
        cls.payload = json.loads(OUT.read_text(encoding="utf-8"))

    def test_complete_affine_transition_matrix(self):
        self.assertEqual(self.payload["status"], "P2726_AFFINE_ORIENTATION_FLOW_NONZERO_BUT_UNSOURCED_NO_STRICT_DYNAMIC_CHIRAL_SOURCE")
        matrix = self.payload["transition_matrix"]
        self.assertEqual(matrix["checked_transition_rows"], 576)
        self.assertEqual(matrix["rule_counts"]["orientation_preserving_translation"], 288)
        self.assertEqual(matrix["rule_counts"]["orientation_flipping_affine_transition"], 288)
        self.assertEqual(matrix["orientation_preserving_delta_values"], [0.0])
        self.assertEqual(matrix["orientation_flipping_delta_values"], [-4.0, 4.0])

    def test_nonzero_flip_is_unsourced_not_closure(self):
        matrix = self.payload["transition_matrix"]
        self.assertEqual(matrix["orientation_preserving_nonzero_count"], 0)
        self.assertEqual(matrix["orientation_flipping_nonzero_count"], 288)
        self.assertTrue(matrix["all_preserving_deltas_zero"])
        self.assertTrue(matrix["all_flipping_deltas_nonzero"])
        acceptance = self.payload["acceptance_matrix"]
        self.assertTrue(acceptance["facts"]["orientation_flip_has_nonzero_signed_jump"])
        self.assertFalse(acceptance["facts"]["strict_orientation_flip_source_law_exported"])
        self.assertFalse(acceptance["accepted_as_strict_dynamic_chiral_source"])

    def test_no_closure_flags_are_exported(self):
        decision = self.payload["decision"]
        self.assertFalse(decision["strict_orientation_flip_source_law_exported"])
        self.assertFalse(decision["nonpremise_orientation_flip_direction_selected"])
        self.assertFalse(decision["p2721_coupling_polarity_selected"])
        self.assertTrue(all(value is False for value in decision["negative_export_flags"].values()))

    def test_docs_written(self):
        self.assertIn("P2726/S1676", MD.read_text(encoding="utf-8"))
        self.assertIn("P2726/S1676", (ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md").read_text(encoding="utf-8"))
        self.assertIn("P2726/S1676", (ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md").read_text(encoding="utf-8"))
        self.assertIn("P2726", (REPO / "AGENTS.md").read_text(encoding="utf-8"))


if __name__ == "__main__":
    unittest.main()
