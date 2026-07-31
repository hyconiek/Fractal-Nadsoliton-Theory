"""Regression and scientific-boundary tests for Programs P365--P378."""

from __future__ import annotations

import csv
import json
import math
from pathlib import Path
import subprocess
import unittest


ROOT = Path(__file__).resolve().parent
RESULTS = json.loads(
    (ROOT / "FIN_Programs_365_378_Results.json").read_text(encoding="utf-8")
)


class Programs365To378Tests(unittest.TestCase):
    def test_summary_contains_all_fourteen_programs(self) -> None:
        with (ROOT / "FIN_Programs_365_378_Summary.csv").open(
            encoding="utf-8"
        ) as handle:
            rows = list(csv.DictReader(handle))
        self.assertEqual(
            [row["program"] for row in rows],
            [f"P{index}" for index in range(365, 379)],
        )

    def test_p365_independent_checker(self) -> None:
        result = RESULTS["P365"]
        self.assertTrue(result["all_checks_pass"])
        self.assertFalse(result["imports_FIN_research_modules"])
        self.assertEqual(result["P351_endpoint_identity_checks"], 12)
        self.assertEqual(result["P352_exact_identity_checks"], 36)
        self.assertIn("[Blocked]", result["status"])

    def test_p366_tightens_both_enclosure_sides(self) -> None:
        result = RESULTS["P366"]
        lower, upper = result["certified_optimum_enclosure"]
        old_lower, old_upper = result["old_enclosure"]
        self.assertGreater(lower, old_lower)
        self.assertLess(upper, old_upper)
        self.assertGreater(result["width_reduction_fraction"], 0.08)
        self.assertEqual(result["certified_atom_count"], 7)
        self.assertEqual(
            result["certified_weight_signs"], [-1, 1, 1, 1, 1, -1, 1]
        )
        self.assertTrue(result["all_krawczyk_components_strictly_inside"])

    def test_external_programs_remain_closed(self) -> None:
        for program in ("P367", "P370", "P372", "P375", "P376", "P378"):
            self.assertFalse(RESULTS[program]["admitted"], program)

    def test_p368_maximal_naturality_audit(self) -> None:
        result = RESULTS["P368"]
        self.assertTrue(result["exact_strict_values_d1_to_d4_pairwise_disjoint"])
        self.assertEqual(result["total_graph_homomorphisms"], 77130)
        self.assertEqual(result["total_isometric_embeddings"], 1664)
        self.assertEqual(result["total_strict_kernel_preserving_maps"], 1664)
        self.assertEqual(result["isometry_kernel_preservation_mismatches"], 0)

    def test_p369_component_identifiability_boundary(self) -> None:
        result = RESULTS["P369"]
        ranks = result["jacobian_audit"]
        self.assertEqual(ranks["complex_transfer"]["effective_rank"], 132)
        self.assertLess(ranks["intensity_only"]["effective_rank"], 132)
        self.assertLess(ranks["conditional_intensity"]["effective_rank"], 132)
        self.assertGreater(
            ranks["conditional_intensity"]["effective_condition_number"], 1e7
        )
        self.assertLess(
            result["worst_tested_setting"]["p95_maximum_vertex_tv"], 0.006
        )
        self.assertIn("synthetic", result["status"])

    def test_p371_analytic_thresholds(self) -> None:
        result = RESULTS["P371"]
        diameter = result["relative_generator_spectral_diameter"]
        thresholds = result["perfect_discrimination_times"]
        for uses in range(1, 5):
            expected = math.pi / (uses * diameter)
            self.assertAlmostEqual(thresholds[str(uses)], expected, places=14)
            self.assertLess(
                abs(
                    result["P343_sampled_threshold_minus_analytic_threshold"][
                        str(uses)
                    ]
                ),
                4e-4,
            )
        self.assertIn("through the first perfect time", result["status"])
        self.assertIn("does not claim", result["boundary"])

    def test_p373_all_five_data_processing_checks(self) -> None:
        result = RESULTS["P373"]
        self.assertEqual(len(result["operational_semantics"]), 5)
        self.assertTrue(
            all(
                violation <= 1e-11
                for violation in result[
                    "maximum_observed_monotonicity_violations"
                ].values()
            )
        )
        self.assertIn("[Conditional]", result["status"])

    def test_p374_outer_extension_is_noncanonical(self) -> None:
        result = RESULTS["P374"]
        self.assertTrue(result["rho_one_half_outer"])
        self.assertGreater(result["certified_rho_threshold_decimal"], 0.66)
        self.assertLess(result["certified_rho_threshold_decimal"], 2 / 3)
        self.assertIn("inserted", result["no_source_result"])

    def test_p377_universality_boundary(self) -> None:
        result = RESULTS["P377"]
        self.assertEqual(result["additivity_failures"], 0)
        self.assertEqual(result["generator_extension_uniqueness_failures"], 0)
        self.assertEqual(result["operational_counterexample_total_variation"], 0.5)
        self.assertIn("[Refuted]", result["status"])

    def test_formal_core_compiles(self) -> None:
        verification = RESULTS["formal_verification"]
        lean = Path(verification["lean_binary"])
        completed = subprocess.run(
            [str(lean), "FIN_Programs_365_378_Formal_Core.lean"],
            cwd=ROOT,
            capture_output=True,
            text=True,
            check=False,
        )
        self.assertEqual(completed.returncode, 0, completed.stderr)

    def test_raw_table_cardinalities(self) -> None:
        expected = {
            "FIN_Programs_365_378_Independent_Checker.csv": 12,
            "FIN_Programs_365_378_Oscillatory_Tightening.csv": 24,
            "FIN_Programs_365_378_Holdout_Gate.csv": 1,
            "FIN_Programs_365_378_Maximal_Naturality.csv": 607,
            "FIN_Programs_365_378_Component_Loss.csv": 18,
            "FIN_Programs_365_378_Hardware_Gate.csv": 1,
            "FIN_Programs_365_378_Subcritical_Comb.csv": 324,
            "FIN_Programs_365_378_Clock_Gate.csv": 1,
            "FIN_Programs_365_378_Resource_Semantics.csv": 320,
            "FIN_Programs_365_378_Outer_Extension.csv": 5,
            "FIN_Programs_365_378_EW_Gate.csv": 1,
            "FIN_Programs_365_378_Standards_Gate.csv": 1,
            "FIN_Programs_365_378_FORC_Universality.csv": 66,
            "FIN_Programs_365_378_Reservoir_Gate.csv": 1,
        }
        for filename, count in expected.items():
            with (ROOT / filename).open(encoding="utf-8") as handle:
                rows = list(csv.DictReader(handle))
            self.assertEqual(len(rows), count, filename)

    def test_all_four_figures_exist(self) -> None:
        expected = {
            "p365_p366_certificate_tightening.png",
            "p368_p369_naturality_identifiability.png",
            "p369_p371_loss_comb.png",
            "p373_p374_semantics_outer.png",
        }
        actual = {
            path.name
            for path in (ROOT / "FIN_Programs_365_378_Figures").glob("*.png")
        }
        self.assertEqual(actual, expected)

    def test_json_has_no_nonstandard_numbers(self) -> None:
        text = (ROOT / "FIN_Programs_365_378_Results.json").read_text(
            encoding="utf-8"
        )
        self.assertNotIn("NaN", text)
        self.assertNotIn("Infinity", text)


if __name__ == "__main__":
    unittest.main()
