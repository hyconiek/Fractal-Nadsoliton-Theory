"""Regression and scientific-boundary tests for Programs P351--P364."""

from __future__ import annotations

import csv
import json
from pathlib import Path
import subprocess
import unittest


ROOT = Path(__file__).resolve().parent
RESULTS = json.loads(
    (ROOT / "FIN_Programs_351_364_Results.json").read_text(encoding="utf-8")
)


class Programs351To364Tests(unittest.TestCase):
    def test_summary_contains_all_fourteen_programs(self) -> None:
        with (ROOT / "FIN_Programs_351_364_Summary.csv").open(
            encoding="utf-8"
        ) as handle:
            rows = list(csv.DictReader(handle))
        self.assertEqual(
            [row["program"] for row in rows],
            [f"P{index}" for index in range(351, 365)],
        )

    def test_p351_krawczyk_certificate_and_enclosure(self) -> None:
        result = RESULTS["P351"]
        self.assertTrue(result["all_krawczyk_components_strictly_inside"])
        self.assertTrue(result["unique_positive_ordered_solution_in_box"])
        lower, upper = result["certified_optimum_enclosure"]
        self.assertLess(lower, upper)
        self.assertLess(upper - lower, 1e-8)
        self.assertGreater(lower, 0.40670632)

    def test_p352_oscillatory_interval_certificate(self) -> None:
        result = RESULTS["P352"]
        self.assertTrue(result["all_twelve_weight_signs_certified"])
        self.assertEqual(result["positive_weight_count"], 8)
        self.assertEqual(result["negative_weight_count"], 4)
        lower, upper = result["certified_optimum_enclosure"]
        self.assertGreater(lower, 0.70727)
        self.assertLess(upper, 0.70737)
        self.assertLess(upper - lower, 1e-4)
        self.assertGreaterEqual(result["safe_dual_bernstein_range"][0], -1.0)
        self.assertLessEqual(result["safe_dual_bernstein_range"][1], 0.0)

    def test_external_gates_remain_closed(self) -> None:
        for program in ("P353", "P356", "P361", "P364"):
            self.assertFalse(RESULTS[program]["admitted"], program)
        self.assertFalse(RESULTS["P356"]["hardware_execution_performed"])
        self.assertFalse(RESULTS["P364"]["hardware_execution_performed"])

    def test_p354_naturality_boundary(self) -> None:
        result = RESULTS["P354"]
        self.assertEqual(result["maximum_isomorphism_naturality_defect"], 0.0)
        self.assertGreater(result["homomorphism_endpoint_defect"], 0.25)
        self.assertIn("[Refuted]", result["status"])

    def test_p355_is_explicitly_synthetic(self) -> None:
        result = RESULTS["P355"]
        self.assertEqual(result["mesh_rotations"], 66)
        self.assertIn("digital twin", result["status"])
        self.assertLess(result["worst_tested_p95_corrected_tv"], 0.002)
        self.assertIn("not a fabricated device", result["boundary"])

    def test_p357_perfect_comb_witnesses(self) -> None:
        result = RESULTS["P357"]
        self.assertEqual(result["uses_tested"], [1, 2, 3, 4])
        self.assertLess(result["maximum_convex_zero_residual"], 1e-12)
        self.assertTrue(result["all_witnesses_have_at_most_three_active_modes"])
        self.assertIn("remains unsolved", result["open_subcritical_problem"])

    def test_p358_clock_design(self) -> None:
        result = RESULTS["P358"]
        self.assertEqual(result["minimum_equal_intervals_for_target"], 7)
        self.assertIn("external", result["external_anchor_required"])

    def test_p359_complete_cancellative_resource_preorder(self) -> None:
        result = RESULTS["P359"]
        self.assertEqual(result["completeness_failures"], 0)
        self.assertEqual(result["cancellation_failures"], 0)
        self.assertEqual(len(result["resource_coordinates"]), 5)
        self.assertEqual(result["ordered_pair_count"], 59049)

    def test_p360_inner_factor_obstruction(self) -> None:
        result = RESULTS["P360"]
        self.assertLess(result["maximum_sampled_inner_factor_modulus_defect"], 1e-12)
        self.assertIn("[Refuted]", result["status"])
        self.assertIn("minimum-phase", result["theorem"])

    def test_p362_metrology_is_conditional(self) -> None:
        result = RESULTS["P362"]
        self.assertEqual(result["exponent_matrix_rank"], 3)
        self.assertTrue(result["benchmark_is_design_only"])
        self.assertFalse(result["external_standards_admitted"])

    def test_p363_axioms_and_boundaries(self) -> None:
        result = RESULTS["P363"]
        self.assertEqual(result["axiom_count"], 5)
        self.assertEqual(result["name"], "FIN Operational Resource Completion (FORC)")
        self.assertIn("not the strict-core selector", result["boundary"])

    def test_formal_core_compiles(self) -> None:
        verification = RESULTS["formal_verification"]
        lean = Path(verification["lean_binary"])
        completed = subprocess.run(
            [str(lean), "FIN_Programs_351_364_Formal_Core.lean"],
            cwd=ROOT,
            capture_output=True,
            text=True,
            check=False,
        )
        self.assertEqual(completed.returncode, 0, completed.stderr)

    def test_raw_table_cardinalities(self) -> None:
        expected = {
            "FIN_Programs_351_364_Envelope_Interval.csv": 6,
            "FIN_Programs_351_364_Oscillatory_Interval.csv": 36,
            "FIN_Programs_351_364_Holdout_Gate.csv": 1,
            "FIN_Programs_351_364_Graph_Naturality.csv": 16,
            "FIN_Programs_351_364_Lossy_Mesh.csv": 16,
            "FIN_Programs_351_364_Hardware_Gate.csv": 1,
            "FIN_Programs_351_364_Adaptive_Comb.csv": 12,
            "FIN_Programs_351_364_Clock_Design.csv": 11,
            "FIN_Programs_351_364_Resource_Category.csv": 128,
            "FIN_Programs_351_364_Analytic_Phase.csv": 27,
            "FIN_Programs_351_364_EW_Gate.csv": 1,
            "FIN_Programs_351_364_Conversion_Metrology.csv": 7,
            "FIN_Programs_351_364_Axiom_Independence.csv": 5,
            "FIN_Programs_351_364_Reservoir_Gate.csv": 1,
        }
        for filename, count in expected.items():
            with (ROOT / filename).open(encoding="utf-8") as handle:
                rows = list(csv.DictReader(handle))
            self.assertEqual(len(rows), count, filename)

    def test_all_four_figures_exist(self) -> None:
        expected = {
            "p351_p352_rigorous_enclosures.png",
            "p354_p355_naturality_loss.png",
            "p357_p358_comb_clock.png",
            "p360_p362_phase_metrology.png",
        }
        actual = {
            path.name
            for path in (ROOT / "FIN_Programs_351_364_Figures").glob("*.png")
        }
        self.assertEqual(actual, expected)

    def test_json_has_no_nonstandard_numbers(self) -> None:
        text = (ROOT / "FIN_Programs_351_364_Results.json").read_text(
            encoding="utf-8"
        )
        self.assertNotIn("NaN", text)
        self.assertNotIn("Infinity", text)


if __name__ == "__main__":
    unittest.main()
