from __future__ import annotations

import json
import subprocess
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent
SCRIPT = ROOT / "p2631_s1581_neural_information_flux_beta_criticality_audit.py"
OUT = ROOT / "generated" / "p2631_s1581_neural_information_flux_beta_criticality_audit.json"
MD = ROOT / "generated" / "p2631_s1581_neural_information_flux_beta_criticality_audit.md"


class P2631NeuralInformationFluxBetaCriticalityAuditTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        subprocess.run([sys.executable, str(SCRIPT)], cwd=ROOT, check=True)
        cls.payload = json.loads(OUT.read_text(encoding="utf-8"))

    def test_content_first_grep_is_present(self) -> None:
        audit = self.payload["semantic_rg_antiduplication_audit"]
        self.assertIn("content-first", audit["mode"])
        self.assertTrue(all(data["count"] > 0 for data in audit["patterns"].values()))

    def test_flux_is_monotone_and_calibratable_not_beta_identity(self) -> None:
        cert = self.payload["flux_conservation_certificate"]
        rows = cert["beta_rows"]
        fluxes = [row["finite_window_flux"] for row in rows]
        self.assertEqual(fluxes, sorted(fluxes, reverse=True))
        self.assertTrue(all(row["flux_derivative"] < 0 for row in rows))
        self.assertIn("choosing C=F(1) makes beta=1 a tautological calibration", cert["consequence"])
        self.assertFalse(cert["exports_beta_equals_1_identity"])

    def test_edge_of_chaos_proxy_does_not_peak_at_beta_one(self) -> None:
        scan = self.payload["edge_of_chaos_scan"]
        self.assertNotEqual(scan["max_entropy_row"]["beta"], 1.0)
        self.assertGreater(scan["max_entropy_row"]["differential_entropy_proxy"], scan["beta_1_row"]["differential_entropy_proxy"])
        self.assertFalse(scan["exports_edge_of_chaos_beta_identity"])

    def test_uv_invariance_blocks_bare_beta_one(self) -> None:
        cert = self.payload["uv_invariance_certificate"]
        self.assertIn("beta * L^eta", cert["invariant"])
        nontrivial = [row for row in cert["rows"] if row["length_unit_scale_a"] != 1.0]
        self.assertTrue(all(not row["bare_beta_equals_1_after_rescaling"] for row in nontrivial))
        self.assertFalse(cert["exports_uv_independent_bare_beta_identity"])

    def test_source_acceptance_and_classification_remain_negative(self) -> None:
        acceptance = self.payload["source_acceptance"]
        self.assertFalse(acceptance["accepts_information_conservation_beta_identity"])
        self.assertIn("flux_conservation_selects_beta_1_without_calibrating_const", acceptance["failed_gates"])
        classification = self.payload["admissible_beta_classification"]
        self.assertIn("admissible", classification["finite_window_norm_stability"]["beta_positive"])
        self.assertTrue(all(value is False for value in self.payload["negative_export_flags"].values()))

    def test_recommendation_and_docs_updated(self) -> None:
        self.assertIn("dimensionless conservation constant", self.payload["recommended_next_honest_step"])
        self.assertIn("P2631/S1581", MD.read_text(encoding="utf-8"))
        self.assertIn("P2631/S1581", (ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md").read_text(encoding="utf-8"))
        self.assertIn("P2631/S1581", (ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md").read_text(encoding="utf-8"))


if __name__ == "__main__":
    unittest.main()
