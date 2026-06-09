from __future__ import annotations

import json
import subprocess
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent
SCRIPT = ROOT / "p2625_s1575_nonlinear_damping_completion_source_classification.py"
OUT = ROOT / "generated" / "p2625_s1575_nonlinear_damping_completion_source_classification.json"
MD = ROOT / "generated" / "p2625_s1575_nonlinear_damping_completion_source_classification.md"


class P2625NonlinearDampingCompletionSourceClassificationTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        subprocess.run([sys.executable, str(SCRIPT)], cwd=ROOT, check=True)
        cls.payload = json.loads(OUT.read_text(encoding="utf-8"))

    def test_content_first_grep_is_present(self) -> None:
        audit = self.payload["semantic_rg_antiduplication_audit"]
        self.assertIn("content-first", audit["mode"])
        self.assertTrue(all(data["count"] > 0 for data in audit["patterns"].values()))

    def test_exact_required_z_beta_isolated(self) -> None:
        self.assertEqual(self.payload["constants"]["beta_tors"], "1/100")
        self.assertEqual(self.payload["constants"]["beta_strict"], "1")
        self.assertEqual(self.payload["constants"]["eta_strict"], "9/5")
        self.assertEqual(self.payload["constants"]["required_Z_beta"], "100")

    def test_candidate_models_reject_shortcuts_and_accept_only_conditional_schema(self) -> None:
        models = {model["name"]: model for model in self.payload["candidate_completion_models"]}
        self.assertFalse(models["legacy_scalar_linear_torsion"]["identity_check"]["identity_for_all_positive_d"])
        self.assertFalse(models["fractal_measure_pushforward_only"]["identity_check"]["identity_for_all_positive_d"])
        self.assertTrue(models["fractal_pushforward_plus_independent_Z_beta"]["identity_check"]["identity_for_all_positive_d"])
        self.assertTrue(models["fractal_pushforward_plus_independent_Z_beta"]["accepted_as_conditional_completion_schema"])
        self.assertEqual(models["fractal_pushforward_plus_independent_Z_beta"]["max_abs_sample_residual"], 0.0)

    def test_source_lattice_requires_positive_beta_source(self) -> None:
        lattice = self.payload["nonlinear_completion_source_lattice"]
        self.assertEqual(lattice["row_count"], 8)
        self.assertEqual(lattice["accepting_row_count"], 1)
        self.assertFalse(lattice["current_unconditional_repair"])
        self.assertEqual(lattice["current_missing_atoms"], ["positive_beta_renormalization_source"])
        self.assertEqual(lattice["minimal_new_required_atom"], "positive_beta_renormalization_source")

    def test_p2620_remains_unrepaired_and_docs_updated(self) -> None:
        self.assertEqual(self.payload["p2620_update"]["accepting_row_count"], 0)
        self.assertTrue(all(value is False for value in self.payload["negative_export_flags"].values()))
        self.assertIn("P2625/S1575", MD.read_text(encoding="utf-8"))
        self.assertIn("P2625/S1575", (ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md").read_text(encoding="utf-8"))
        self.assertIn("P2625/S1575", (ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md").read_text(encoding="utf-8"))


if __name__ == "__main__":
    unittest.main()
