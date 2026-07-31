#!/usr/bin/env python3
"""Regression tests for FIN Research Programs P379--P393."""

from __future__ import annotations

import csv
import json
import math
from pathlib import Path
import subprocess
import unittest


ROOT = Path(__file__).resolve().parent
RESULTS = json.loads(
    (ROOT / "FIN_Programs_379_393_Results.json").read_text(encoding="utf-8")
)
LEAN = (
    ROOT
    / ".elan"
    / "toolchains"
    / "leanprover--lean4---v4.28.0"
    / "bin"
    / "lean"
)


def csv_rows(filename: str) -> list[dict[str, str]]:
    with (ROOT / filename).open(encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


class Programs379To393Tests(unittest.TestCase):
    def test_summary_contains_all_fifteen_programs(self) -> None:
        rows = csv_rows("FIN_Programs_379_393_Summary.csv")
        self.assertEqual([row["program"] for row in rows], [
            f"P{program}" for program in range(379, 394)
        ])

    def test_p379_reflects_thirty_five_predicates(self) -> None:
        result = RESULTS["P379"]
        self.assertEqual(result["kernel_checked_predicate_count"], 35)
        self.assertEqual(result["lean_returncode"], 0)
        self.assertEqual(
            len(csv_rows("FIN_Programs_379_393_Arithmetic_Reflection.csv")),
            35,
        )

    def test_both_lean_sources_compile(self) -> None:
        for source, timeout in (
            ("FIN_Programs_379_393_Formal_Core.lean", 60),
            ("FIN_Programs_379_393_Arithmetic_Reflection.lean", 180),
        ):
            completed = subprocess.run(
                [str(LEAN), source],
                cwd=ROOT,
                capture_output=True,
                text=True,
                check=False,
                timeout=timeout,
            )
            self.assertEqual(completed.returncode, 0, completed.stderr)

    def test_p380_gap_is_below_one_e_minus_seven(self) -> None:
        result = RESULTS["P380"]
        self.assertTrue(result["gap_below_1e_7"])
        self.assertLess(result["certified_gap"], 1e-7)
        lower, upper = result["certified_enclosure"]
        old_lower, old_upper = result["previous_P366_enclosure"]
        self.assertGreater(lower, old_lower)
        self.assertLessEqual(upper, old_upper)
        self.assertLessEqual(lower, upper)

    def test_p381_global_injectivity_certificate(self) -> None:
        result = RESULTS["P381"]
        self.assertEqual(result["finite_double_collision_count"], 0)
        self.assertIn("Lindemann-Weierstrass", result["theorem"])
        self.assertEqual(result["domain"], "d in N_0")

    def test_p382_is_one_target_defined_atom_only(self) -> None:
        result = RESULTS["P382"]
        rows = csv_rows("FIN_Programs_379_393_Damping_Naturality.csv")
        self.assertFalse(result["source_independent"])
        self.assertGreater(
            result["C12_full_kernel_relative_l2_residual_after_damping_only"],
            0.9,
        )
        self.assertLess(
            max(abs(float(row["attenuation_identity_defect"])) for row in rows),
            1e-15,
        )

    def test_p383_coordinate_gauge_boundary(self) -> None:
        result = RESULTS["P383"]
        self.assertEqual(result["component_count"], 66)
        self.assertLess(result["maximum_combined_lattice_defect"], 1e-12)
        self.assertGreater(result["minimum_single_pi_defect"], 0.1)
        self.assertGreater(result["minimum_single_loss_shift_defect"], 1e-4)
        self.assertIn("[Open]", result["status"])

    def test_external_programs_remain_closed(self) -> None:
        for program in ("P384", "P390", "P391", "P392", "P393"):
            self.assertFalse(RESULTS[program]["admitted"])
            self.assertIn("[Blocked by external evidence]", RESULTS[program]["status"])

    def test_p385_all_feasible_lowers_respect_upper(self) -> None:
        rows = csv_rows("FIN_Programs_379_393_Lossy_Comb.csv")
        for row in rows:
            self.assertLessEqual(
                float(row["best_feasible_lower"]),
                float(row["adaptive_hybrid_upper"]) + 2e-12,
            )
        ideal = [
            row for row in rows
            if row["survival_probability"] == "1.0"
            and row["coherence_factor"] == "1.0"
            and row["fraction_of_ideal_threshold"] == "1.0"
        ]
        self.assertEqual(len(ideal), 4)
        for row in ideal:
            self.assertTrue(
                math.isclose(float(row["ghz_lower"]), 1.0, abs_tol=1e-12)
            )

    def test_p386_clock_tubes_are_ordered(self) -> None:
        rows = csv_rows("FIN_Programs_379_393_Clock_Tube.csv")
        for row in rows:
            self.assertLessEqual(float(row["tau_lower"]), float(row["tau_upper"]))
            self.assertLessEqual(
                float(row["discrimination_lower"]),
                float(row["discrimination_upper"]) + 1e-15,
            )
        self.assertGreater(RESULTS["P386"]["maximum_reported_discrimination_band_width"], 0)

    def test_p387_realization_is_a_probability_model(self) -> None:
        rows = csv_rows("FIN_Programs_379_393_Jordan_Realization.csv")
        probability = sum(float(row["sampling_probability"]) for row in rows)
        self.assertTrue(math.isclose(probability, 1.0, abs_tol=2e-15))
        self.assertLess(RESULTS["P387"]["maximum_center_moment_residual"], 2e-15)
        self.assertTrue(
            math.isclose(
                RESULTS["P387"]["negative_mass"],
                RESULTS["P380"]["certified_enclosure"][1],
                abs_tol=2e-15,
            )
        )

    def test_p388_refutes_all_twenty_cross_implications(self) -> None:
        result = RESULTS["P388"]
        self.assertEqual(result["identity_conversions"], 5)
        self.assertEqual(
            result["directed_cross_conversion_implications_refuted"], 20
        )
        self.assertEqual(
            sum(
                row["counterexample"] == "True"
                for row in csv_rows("FIN_Programs_379_393_Conversion_Matrix.csv")
            ),
            20,
        )

    def test_p389_dilation_semigroup(self) -> None:
        result = RESULTS["P389"]
        self.assertLess(result["maximum_truncated_numerical_defect"], 1e-15)
        self.assertIn("rho=rho/2", result["no_go"])
        self.assertGreater(result["certified_outer_domain_upper"], 0.66)

    def test_all_four_figures_exist(self) -> None:
        expected = {
            "p379_p380_reflection_contact.png",
            "p381_p383_injectivity_gauge.png",
            "p385_p386_loss_clock.png",
            "p387_p389_realization_scale.png",
        }
        actual = {
            path.name for path in (ROOT / "FIN_Programs_379_393_Figures").glob("*.png")
        }
        self.assertEqual(actual, expected)
        for filename in expected:
            self.assertGreater(
                (ROOT / "FIN_Programs_379_393_Figures" / filename).stat().st_size,
                10_000,
            )

    def test_raw_table_cardinalities(self) -> None:
        expected = {
            "FIN_Programs_379_393_Arithmetic_Reflection.csv": 35,
            "FIN_Programs_379_393_Oscillatory_Contact.csv": 12,
            "FIN_Programs_379_393_Distance_Injectivity.csv": 65,
            "FIN_Programs_379_393_Damping_Naturality.csv": 7,
            "FIN_Programs_379_393_Component_Gauge.csv": 66,
            "FIN_Programs_379_393_Photonic_Gate.csv": 1,
            "FIN_Programs_379_393_Lossy_Comb.csv": 192,
            "FIN_Programs_379_393_Clock_Tube.csv": 16,
            "FIN_Programs_379_393_Jordan_Realization.csv": 7,
            "FIN_Programs_379_393_Conversion_Matrix.csv": 25,
            "FIN_Programs_379_393_Outer_Scale.csv": 12,
            "FIN_Programs_379_393_Holdout_Gate.csv": 1,
            "FIN_Programs_379_393_Standards_Gate.csv": 1,
            "FIN_Programs_379_393_Reservoir_Gate.csv": 1,
            "FIN_Programs_379_393_EW_Gate.csv": 1,
        }
        for filename, count in expected.items():
            self.assertEqual(len(csv_rows(filename)), count, filename)

    def test_json_has_no_nonstandard_numbers(self) -> None:
        text = (ROOT / "FIN_Programs_379_393_Results.json").read_text(
            encoding="utf-8"
        )
        self.assertNotIn("NaN", text)
        self.assertNotIn("Infinity", text)


if __name__ == "__main__":
    unittest.main()
