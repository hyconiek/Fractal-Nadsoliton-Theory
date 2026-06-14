from __future__ import annotations

import json
import subprocess
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
SCRIPT = ROOT / "p2723_s1673_strict_chiral_time_orientation_source_law_candidate_matrix.py"
OUT = ROOT / "generated" / "p2723_s1673_strict_chiral_time_orientation_source_law_candidate_matrix.json"
MD = ROOT / "generated" / "p2723_s1673_strict_chiral_time_orientation_source_law_candidate_matrix.md"


class P2723StrictChiralTimeOrientationSourceLawCandidateMatrixTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        subprocess.run([sys.executable, str(SCRIPT)], cwd=ROOT, check=True)
        cls.payload = json.loads(OUT.read_text(encoding="utf-8"))

    def test_no_candidate_exports_strict_chiral_time_source(self):
        self.assertEqual(self.payload["status"], "P2723_STRICT_CHIRAL_TIME_ORIENTATION_SOURCE_LAW_MATRIX_NO_ACCEPTED_SOURCE")
        matrix = self.payload["source_law_matrix"]
        self.assertEqual(matrix["accepted_candidate_count"], 0)
        self.assertFalse(matrix["strict_chiral_time_orientation_source_law_exported"])
        self.assertEqual(len(matrix["candidate_rows"]), 4)
        self.assertTrue(all(not row["accepted_as_strict_chiral_time_source"] for row in matrix["candidate_rows"]))

    def test_expected_blockers_are_recorded(self):
        rows = {row["name"]: row for row in self.payload["source_law_matrix"]["candidate_rows"]}
        self.assertIn("premise/convention", rows["external_time_arrow_or_clockwise_premise"]["blocker"])
        self.assertIn("time-series/EOM", rows["chiral_bispectrum_temporal_phase_velocity"]["blocker"])
        self.assertIn("lambda sign paired", rows["inversion_odd_character_lambda_sign_law"]["blocker"])
        self.assertIn("no non-premise source", rows["boundary_cocycle_orientation_flow"]["blocker"])

    def test_no_closure_flags_are_exported(self):
        decision = self.payload["decision"]
        self.assertFalse(decision["strict_chiral_time_orientation_source_law_exported"])
        self.assertFalse(decision["nonzero_signed_time_orientation_value_exported"])
        self.assertFalse(decision["canonical_coupling_polarity_selected"])
        self.assertFalse(decision["strict_mechanism_fixing_lambda_exported"])
        self.assertTrue(all(value is False for value in decision["negative_export_flags"].values()))

    def test_docs_written(self):
        self.assertIn("P2723/S1673", MD.read_text(encoding="utf-8"))
        self.assertIn("P2723/S1673", (ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md").read_text(encoding="utf-8"))
        self.assertIn("P2723/S1673", (ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md").read_text(encoding="utf-8"))
        self.assertIn("P2723", (REPO / "AGENTS.md").read_text(encoding="utf-8"))


if __name__ == "__main__":
    unittest.main()
