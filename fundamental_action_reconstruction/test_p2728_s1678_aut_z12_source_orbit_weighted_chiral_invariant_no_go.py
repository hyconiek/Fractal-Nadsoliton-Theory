from __future__ import annotations

import json
import subprocess
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
SCRIPT = ROOT / "p2728_s1678_aut_z12_source_orbit_weighted_chiral_invariant_no_go.py"
OUT = ROOT / "generated" / "p2728_s1678_aut_z12_source_orbit_weighted_chiral_invariant_no_go.json"
MD = ROOT / "generated" / "p2728_s1678_aut_z12_source_orbit_weighted_chiral_invariant_no_go.md"


class P2728AutZ12SourceOrbitWeightedChiralInvariantNoGoTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        subprocess.run([sys.executable, str(SCRIPT)], cwd=ROOT, check=True)
        cls.payload = json.loads(OUT.read_text(encoding="utf-8"))

    def test_exhausts_aut_orbit_weightings(self):
        self.assertEqual(self.payload["status"], "P2728_AUT_Z12_SOURCE_ORBIT_WEIGHTED_CHIRAL_INVARIANT_NO_GO")
        audit = self.payload["source_orbit_weight_audit"]
        self.assertEqual(audit["orbit_count"], 6)
        self.assertEqual(audit["weight_values"], [-1, 0, 1])
        self.assertEqual(audit["weighting_count"], 729)
        self.assertEqual(audit["nonzero_global_signed_total_count"], 0)
        self.assertTrue(audit["all_global_signed_totals_zero"])

    def test_no_single_polarity_weighting_exists(self):
        audit = self.payload["source_orbit_weight_audit"]
        self.assertEqual(audit["single_polarity_weighting_count"], 0)
        self.assertTrue(audit["all_nonzero_weightings_have_paired_row_polarities"])
        acceptance = self.payload["acceptance_matrix"]
        self.assertFalse(acceptance["facts"]["single_polarity_weighting_exported"])
        self.assertFalse(acceptance["accepted_as_strict_source_dependent_invariant"])

    def test_no_closure_flags_are_exported(self):
        decision = self.payload["decision"]
        self.assertFalse(decision["aut_orbit_weighted_source_invariant_exported"])
        self.assertFalse(decision["nonzero_global_signed_value_exported"])
        self.assertFalse(decision["p2721_coupling_polarity_selected"])
        self.assertTrue(all(value is False for value in decision["negative_export_flags"].values()))

    def test_docs_written(self):
        self.assertIn("P2728/S1678", MD.read_text(encoding="utf-8"))
        self.assertIn("P2728/S1678", (ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md").read_text(encoding="utf-8"))
        self.assertIn("P2728/S1678", (ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md").read_text(encoding="utf-8"))
        self.assertIn("P2728", (REPO / "AGENTS.md").read_text(encoding="utf-8"))


if __name__ == "__main__":
    unittest.main()
