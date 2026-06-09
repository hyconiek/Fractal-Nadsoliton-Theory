from __future__ import annotations

import json
import subprocess
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent
SCRIPT = ROOT / "p2626_s1576_micro_zbeta_source_nonpromotion_audit.py"
OUT = ROOT / "generated" / "p2626_s1576_micro_zbeta_source_nonpromotion_audit.json"
MD = ROOT / "generated" / "p2626_s1576_micro_zbeta_source_nonpromotion_audit.md"


class P2626MicroZbetaSourceNonpromotionAuditTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        subprocess.run([sys.executable, str(SCRIPT)], cwd=ROOT, check=True)
        cls.payload = json.loads(OUT.read_text(encoding="utf-8"))

    def test_content_first_grep_is_present(self) -> None:
        audit = self.payload["semantic_rg_antiduplication_audit"]
        self.assertIn("content-first", audit["mode"])
        self.assertTrue(all(data["count"] > 0 for data in audit["patterns"].values()))

    def test_qw2064_metrics_are_extracted(self) -> None:
        metrics = self.payload["qw2064_metrics"]
        self.assertEqual(metrics["target_z_beta"], 100.0)
        self.assertGreater(metrics["micro_global_median_z_beta"], 0.0)
        self.assertTrue(metrics["target_inside_q25_q75"])
        self.assertTrue(metrics["wide_ci_warning"])
        self.assertTrue(metrics["target_definition_depends_on_frozen_kernel"])
        self.assertFalse(metrics["exact_target_match"])

    def test_positive_beta_source_is_not_promoted(self) -> None:
        table = self.payload["positive_beta_source_acceptance"]
        self.assertFalse(table["accepts_positive_beta_renormalization_source"])
        self.assertIn("target_independent_of_selected_kernel", table["failed_gates"])
        self.assertIn("exact_or_bridge_tolerance_theorem", table["failed_gates"])
        self.assertIn("narrow_dispersion_no_wide_ci_warning", table["failed_gates"])

    def test_p2625_and_ltotal_remain_closed(self) -> None:
        self.assertFalse(self.payload["p2625_update"]["positive_beta_renormalization_source_after_p2626"])
        self.assertFalse(self.payload["p2625_update"]["nonlinear_damping_completion_source_after_p2626"])
        self.assertTrue(all(value is False for value in self.payload["negative_export_flags"].values()))
        self.assertIn("P2626/S1576", MD.read_text(encoding="utf-8"))
        self.assertIn("P2626/S1576", (ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md").read_text(encoding="utf-8"))
        self.assertIn("P2626/S1576", (ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md").read_text(encoding="utf-8"))


if __name__ == "__main__":
    unittest.main()
