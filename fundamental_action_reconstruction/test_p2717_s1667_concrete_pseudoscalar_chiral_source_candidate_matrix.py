from __future__ import annotations

import json
import subprocess
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
SCRIPT = ROOT / "p2717_s1667_concrete_pseudoscalar_chiral_source_candidate_matrix.py"
OUT = ROOT / "generated" / "p2717_s1667_concrete_pseudoscalar_chiral_source_candidate_matrix.json"
MD = ROOT / "generated" / "p2717_s1667_concrete_pseudoscalar_chiral_source_candidate_matrix.md"


class P2717ConcretePseudoscalarChiralSourceCandidateMatrixTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        subprocess.run([sys.executable, str(SCRIPT)], cwd=ROOT, check=True)
        cls.payload = json.loads(OUT.read_text(encoding="utf-8"))

    def test_candidate_matrix_has_no_accepted_source(self):
        self.assertEqual(self.payload["status"], "P2717_CONCRETE_PSEUDOSCALAR_CHIRAL_SOURCE_MATRIX_NO_ACCEPTED_SOURCE")
        self.assertEqual(self.payload["accepted_candidate_count"], 0)
        candidates = {row["candidate_id"] for row in self.payload["candidate_matrix"]}
        self.assertEqual(candidates, {
            "levi_civita_volume_orientation_density",
            "pontryagin_or_chiral_anomaly_density",
            "eta_or_spectral_asymmetry_index",
            "oriented_z12_cycle_cup_product",
        })
        self.assertTrue(all(not row["accepted_as_strict_source"] for row in self.payload["candidate_matrix"]))
        self.assertTrue(all(not row["fixes_lambda_now"] for row in self.payload["candidate_matrix"]))

    def test_each_candidate_has_a_strict_blocker_and_no_closure_flags(self):
        for row in self.payload["candidate_matrix"]:
            self.assertTrue(row["inversion_odd_representation"])
            self.assertTrue(row["missing_criteria"])
        z12 = next(row for row in self.payload["candidate_matrix"] if row["candidate_id"] == "oriented_z12_cycle_cup_product")
        self.assertIn("nonzero_signed_value_exported", z12["missing_criteria"])
        self.assertIn("sign_not_orientation_convention", z12["missing_criteria"])
        decision = self.payload["decision"]
        self.assertFalse(decision["concrete_pseudoscalar_source_accepted"])
        self.assertFalse(decision["strict_mechanism_fixing_lambda_exported"])
        self.assertTrue(all(value is False for value in decision["negative_export_flags"].values()))

    def test_docs_written(self):
        self.assertIn("P2717/S1667", MD.read_text(encoding="utf-8"))
        self.assertIn("P2717/S1667", (ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md").read_text(encoding="utf-8"))
        self.assertIn("P2717/S1667", (ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md").read_text(encoding="utf-8"))
        self.assertIn("P2717", (REPO / "AGENTS.md").read_text(encoding="utf-8"))


if __name__ == "__main__":
    unittest.main()
