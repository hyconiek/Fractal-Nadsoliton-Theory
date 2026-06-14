from __future__ import annotations

import json
import subprocess
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
SCRIPT = ROOT / "p2727_s1677_orientation_transition_law_equivariance_and_polarity_no_go.py"
OUT = ROOT / "generated" / "p2727_s1677_orientation_transition_law_equivariance_and_polarity_no_go.json"
MD = ROOT / "generated" / "p2727_s1677_orientation_transition_law_equivariance_and_polarity_no_go.md"


class P2727OrientationTransitionLawEquivarianceAndPolarityNoGoTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        subprocess.run([sys.executable, str(SCRIPT)], cwd=ROOT, check=True)
        cls.payload = json.loads(OUT.read_text(encoding="utf-8"))

    def test_exhausts_source_independent_orientation_laws(self):
        self.assertEqual(self.payload["status"], "P2727_ORIENTATION_TRANSITION_LAW_EQUIVARIANCE_POLARITY_NO_GO")
        audit = self.payload["orientation_law_audit"]
        self.assertEqual(audit["law_count"], 4)
        self.assertEqual(audit["velocity_count"], 12)
        self.assertEqual(audit["checked_transition_rows"], 1152)
        self.assertEqual(audit["equivariant_law_names"], ["preserve_orientation", "flip_orientation"])

    def test_polarity_selectors_are_premise_laws(self):
        audit = self.payload["orientation_law_audit"]
        self.assertEqual(audit["equivariant_polarity_selecting_law_count"], 0)
        self.assertEqual(
            audit["premise_polarity_selecting_law_names"],
            ["collapse_to_plus_orientation", "collapse_to_minus_orientation"],
        )
        summaries = {row["law"]: row for row in audit["law_summaries"]}
        self.assertEqual(summaries["preserve_orientation"]["delta_values"], [0.0])
        self.assertTrue(summaries["flip_orientation"]["balanced_nonzero_polarities"])

    def test_no_closure_flags_are_exported(self):
        decision = self.payload["decision"]
        self.assertFalse(decision["equivariant_polarity_selecting_orientation_law_exported"])
        self.assertFalse(decision["nonpremise_orientation_branch_selected"])
        self.assertFalse(decision["p2721_coupling_polarity_selected"])
        self.assertTrue(all(value is False for value in decision["negative_export_flags"].values()))

    def test_docs_written(self):
        self.assertIn("P2727/S1677", MD.read_text(encoding="utf-8"))
        self.assertIn("P2727/S1677", (ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md").read_text(encoding="utf-8"))
        self.assertIn("P2727/S1677", (ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md").read_text(encoding="utf-8"))
        self.assertIn("P2727", (REPO / "AGENTS.md").read_text(encoding="utf-8"))


if __name__ == "__main__":
    unittest.main()
