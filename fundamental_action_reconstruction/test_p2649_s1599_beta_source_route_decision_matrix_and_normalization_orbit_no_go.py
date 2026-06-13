from __future__ import annotations

import json
import subprocess
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent
SCRIPT = ROOT / "p2649_s1599_beta_source_route_decision_matrix_and_normalization_orbit_no_go.py"
OUT = ROOT / "generated" / "p2649_s1599_beta_source_route_decision_matrix_and_normalization_orbit_no_go.json"
MD = ROOT / "generated" / "p2649_s1599_beta_source_route_decision_matrix_and_normalization_orbit_no_go.md"


class P2649BetaSourceRouteDecisionMatrixTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        subprocess.run([sys.executable, str(SCRIPT)], cwd=ROOT, check=True)
        cls.payload = json.loads(OUT.read_text(encoding="utf-8"))

    def test_content_first_grep_and_upstream_atoms(self) -> None:
        audit = self.payload["semantic_rg_antiduplication_audit"]
        self.assertIn("content-first", audit["mode"])
        self.assertIn("beta_source_identity_content", audit["patterns"])
        self.assertIn("normalization_orbit_content", audit["patterns"])
        self.assertIn("micro_zbeta_mismatch_content", audit["patterns"])
        self.assertIn("empirical_not_source_content", audit["patterns"])
        self.assertTrue(all(data["count"] > 0 for data in audit["patterns"].values()))
        atoms = self.payload["upstream_atoms"]
        self.assertTrue(atoms["holdout_harness_ready"])
        self.assertTrue(atoms["statistical_margin_rule_ready"])
        self.assertFalse(atoms["holdout_real_blind_data_executed"])
        self.assertFalse(atoms["empirical_holdout_can_replace_source_theorem"])

    def test_normalization_orbit_gauges_all_positive_betas_to_one(self) -> None:
        orbit = self.payload["normalization_orbit_witness"]
        self.assertFalse(orbit["exports_beta_source"])
        self.assertTrue(orbit["all_positive_betas_gauge_to_one"])
        self.assertIn("beta'=1", orbit["theorem"])
        for row in orbit["rows"]:
            self.assertAlmostEqual(row["beta_prime_after_rescaling"], 1.0, places=12)
            self.assertTrue(row["bare_beta_one_reached_by_rescaling"])

    def test_tail_ratio_inversion_recovers_beta_but_does_not_source_it(self) -> None:
        witness = self.payload["tail_ratio_inversion_witness"]
        self.assertFalse(witness["exports_beta_source"])
        self.assertTrue(witness["all_strict_ratios_recover_beta_one"])
        for row in witness["rows"]:
            self.assertAlmostEqual(row["beta_recovered_from_strict_ratio"], 1.0, places=12)
            self.assertAlmostEqual(row["beta_recovered_from_legacy_ratio"], 0.01, places=12)

    def test_no_beta_source_route_passes_and_docs_updated(self) -> None:
        matrix = self.payload["beta_source_route_matrix"]
        self.assertEqual(len(matrix), 5)
        self.assertFalse(any(row["passes_as_target_independent_beta_source"] for row in matrix))
        by_route = {row["route"]: row for row in matrix}
        self.assertIn("normalization_gauge_fixed", by_route["normalization_orbit_beta_one"]["missing_required_atoms"])
        self.assertIn("target_independent_conservation_constant_exported", by_route["information_flux_conservation_beta_one"]["missing_required_atoms"])
        self.assertIn("positive_zbeta_source_exported", by_route["micro_zbeta_bridge_source"]["missing_required_atoms"])
        self.assertIn("empirical_holdout_can_replace_source_theorem", by_route["blind_empirical_compression_holdout"]["missing_required_atoms"])
        decision = self.payload["closure_decision"]
        self.assertEqual(decision["passing_beta_source_routes"], [])
        self.assertFalse(decision["beta_source_exported_now"])
        self.assertFalse(decision["full_kernel_now"])
        self.assertFalse(decision["toe_closure_now"])
        self.assertTrue(all(value is False for value in self.payload["negative_export_flags"].values()))
        self.assertIn("P2649/S1599", MD.read_text(encoding="utf-8"))
        self.assertIn("P2649/S1599", (ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md").read_text(encoding="utf-8"))
        self.assertIn("P2649/S1599", (ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md").read_text(encoding="utf-8"))


if __name__ == "__main__":
    unittest.main()
