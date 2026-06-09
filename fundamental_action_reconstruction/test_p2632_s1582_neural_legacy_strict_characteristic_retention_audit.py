from __future__ import annotations

import json
import subprocess
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent
SCRIPT = ROOT / "p2632_s1582_neural_legacy_strict_characteristic_retention_audit.py"
OUT = ROOT / "generated" / "p2632_s1582_neural_legacy_strict_characteristic_retention_audit.json"
MD = ROOT / "generated" / "p2632_s1582_neural_legacy_strict_characteristic_retention_audit.md"


class P2632NeuralLegacyStrictCharacteristicRetentionAuditTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        subprocess.run([sys.executable, str(SCRIPT)], cwd=ROOT, check=True)
        cls.payload = json.loads(OUT.read_text(encoding="utf-8"))

    def test_content_first_grep_is_present(self) -> None:
        audit = self.payload["semantic_rg_antiduplication_audit"]
        self.assertIn("content-first", audit["mode"])
        self.assertTrue(all(data["count"] > 0 for data in audit["patterns"].values()))

    def test_phase_is_retained_but_denominator_is_strict_addition(self) -> None:
        geometry = self.payload["finite_sample_geometry"]
        self.assertEqual(geometry["phase_linf_difference"], 0.0)
        self.assertLess(geometry["strict_attention_at_d10_over_legacy"], 0.02)
        self.assertGreater(geometry["strict_slope_at_d10_minus_legacy_slope"], 1.0)
        self.assertIn("not a small perturbation", geometry["interpretation"])

    def test_characteristic_matrix_does_not_claim_full_preservation(self) -> None:
        matrix = self.payload["characteristic_retention_matrix"]
        self.assertEqual(len(matrix), 6)
        self.assertTrue(any(row["strict_retention_verdict"] == "enriched_strict_successor_not_exact_retention" for row in matrix))
        self.assertTrue(any(row["strict_retention_verdict"] == "not_preserved_as_bare_parameter" for row in matrix))
        self.assertTrue(all(not row["preserved_as_final"] for row in matrix))

    def test_professorial_readout_marks_working_incomplete_kernel(self) -> None:
        readout = self.payload["professorial_readout"]
        self.assertFalse(readout["does_strict_preserve_all_legacy_characteristics"])
        self.assertEqual(readout["status"], "STRICT_KERNEL_IS_ROBUST_WORKING_SUCCESSOR_BUT_INCOMPLETE_FOR_TOE_FINALITY")
        self.assertTrue(any("beta_tors=0.01 is not the strict beta=1" in item for item in readout["why_not"]))

    def test_no_finality_or_role_transfer_export(self) -> None:
        acceptance = self.payload["source_acceptance"]
        self.assertFalse(acceptance["accepts_strict_kernel_as_complete_toe_kernel"])
        self.assertIn("beta_tors_to_beta_1_bridge_source_closed", acceptance["failed_gates"])
        self.assertTrue(all(value is False for value in self.payload["negative_export_flags"].values()))

    def test_recommendation_and_docs_updated(self) -> None:
        self.assertIn("characteristic-by-characteristic completion certificate", self.payload["recommended_next_honest_step"])
        self.assertIn("P2632/S1582", MD.read_text(encoding="utf-8"))
        self.assertIn("P2632/S1582", (ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md").read_text(encoding="utf-8"))
        self.assertIn("P2632/S1582", (ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md").read_text(encoding="utf-8"))


if __name__ == "__main__":
    unittest.main()
