#!/usr/bin/env python3
"""Regression tests for the local FIN P428--P430 checkpoint."""

from __future__ import annotations

import csv
import json
import math
from pathlib import Path
import subprocess
import unittest


ROOT = Path(__file__).resolve().parent
RESULTS = json.loads(
    (ROOT / "FIN_Programs_428_430_Results.json").read_text(encoding="utf-8")
)
LEAN = ROOT / ".elan/toolchains/leanprover--lean4---v4.28.0/bin/lean"


def table(name: str) -> list[dict[str, str]]:
    with (ROOT / name).open(encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


class Programs428To430Tests(unittest.TestCase):
    def test_p428_type_repair_and_lean_reduction(self) -> None:
        p428 = RESULTS["P428"]
        self.assertEqual(p428["lean_returncode"], 0)
        self.assertEqual(p428["angle_count"], 12)
        self.assertFalse(p428["old_interface_type_adequate_for_real_cosine"])
        self.assertFalse(p428["mathlib_present"])
        self.assertLess(p428["maximum_rational_width"], 1e-30)
        completed = subprocess.run(
            [str(LEAN), "FIN_Programs_428_430_P428_Cosine_Reduction.lean"],
            cwd=ROOT,
            capture_output=True,
            text=True,
            timeout=120,
            check=False,
        )
        self.assertEqual(completed.returncode, 0, completed.stderr)

    def test_p428_all_rational_side_conditions(self) -> None:
        rows = table("FIN_Programs_428_430_P428_Cosine.csv")
        self.assertEqual(len(rows), 12)
        self.assertTrue(all(row["nonnegative"] == "True" for row in rows))
        self.assertTrue(all(row["square_below_12"] == "True" for row in rows))
        self.assertTrue(all(float(row["width"]) > 0 for row in rows))

    def test_p429_exact_inclusion(self) -> None:
        p429 = RESULTS["P429"]
        self.assertTrue(p429["all_25_variables_strictly_included"])
        self.assertEqual(p429["successful_relative_radius_exponent"], 20)
        self.assertEqual(p429["weight_signs"], [-1, 1, 1, 1, 1, -1, 1])
        self.assertTrue(p429["node_order_certified"])
        self.assertLess(p429["weighted_infinity_contraction_upper"], 1.0)
        self.assertLess(float(p429["refinement"]["maximum_residual"]), 1e-140)
        self.assertGreater(p429["refinement"]["double_condition_number"], 1e16)

        rows = table("FIN_Programs_428_430_P429_Krawczyk.csv")
        radius = [row for row in rows if row["row_type"] == "radius_attempt"]
        variables = [row for row in rows if row["row_type"] == "certified_variable"]
        self.assertEqual(len(radius), 1)
        self.assertEqual(radius[0]["all_inside"], "True")
        self.assertEqual(radius[0]["strict_contraction"], "True")
        self.assertEqual(int(radius[0]["inside_variables"]), 25)
        self.assertEqual(len(variables), 25)
        self.assertTrue(all(row["strictly_inside"] == "True" for row in variables))

    def test_p430_global_dual_feasibility(self) -> None:
        p430 = RESULTS["P430"]
        self.assertTrue(p430["global_feasibility_certified"])
        self.assertEqual(p430["stationary_contacts_certified"], 6)
        self.assertTrue(p430["endpoint_certified"])
        self.assertEqual(
            p430["noncontact_gaps_certified"], p430["noncontact_gaps_total"]
        )
        rows = table("FIN_Programs_428_430_P430_Global_Dual.csv")
        self.assertEqual(len(rows), 14)
        self.assertTrue(all(row["safe"] == "True" for row in rows))
        gap_rows = [row for row in rows if row["row_type"] == "noncontact_bernstein_cell"]
        self.assertEqual(len(gap_rows), 7)
        self.assertTrue(all(float(row["lower"]) > -1 for row in gap_rows))
        self.assertTrue(all(float(row["upper"]) < 0 for row in gap_rows))

    def test_objective_and_preconditioner_boundaries(self) -> None:
        p429 = RESULTS["P429"]
        self.assertTrue(
            math.isclose(p429["objective_midpoint"], 0.7073534677231137, abs_tol=1e-15)
        )
        self.assertLessEqual(
            p429["objective_interval_lower"], p429["objective_midpoint"]
        )
        self.assertGreaterEqual(
            p429["objective_interval_upper"], p429["objective_midpoint"]
        )
        self.assertLess(p429["objective_interval_width"], 1e-18)
        self.assertIn("rounded matrix", p429["preconditioner"])
        self.assertIn("exactly", p429["preconditioner"])

    def test_legacy_star_context_is_guarded(self) -> None:
        context = RESULTS["legacy_star_coupling_context"]
        self.assertTrue(all(context["source_checks"].values()))
        self.assertEqual(context["path_sum_exponent_claim"], 1.0)
        self.assertEqual(
            context["path_sum_exponent_required_for_inverse_linear_tail"], -2.6
        )
        self.assertEqual(context["exact_real_zero_set"], "d=4/3+4n")
        self.assertFalse(context["selector_discharged"])
        self.assertFalse(context["physical_role_transfer_exported"])
        alleged = context["alleged_integer_node_cosine_values"]
        self.assertTrue(all(abs(float(value)) > 0.2 for value in alleged.values()))

    def test_global_boundaries(self) -> None:
        metadata = RESULTS["metadata"]
        self.assertTrue(metadata["local_only"])
        self.assertFalse(metadata["external_physical_evidence"])
        self.assertTrue(metadata["kernel_split_preserved"])
        self.assertFalse(metadata["selector_discharged"])
        self.assertFalse(metadata["dimensional_source_exported"])

    def test_json_is_strict_and_summary_complete(self) -> None:
        raw = (ROOT / "FIN_Programs_428_430_Results.json").read_text(encoding="utf-8")
        self.assertNotIn("NaN", raw)
        self.assertNotIn("Infinity", raw)
        summary = table("FIN_Programs_428_430_Summary.csv")
        self.assertEqual([row["program"] for row in summary], ["P428", "P429", "P430"])
        figure = ROOT / "FIN_Programs_428_430_Figures" / "p428_p430_exact_closure_and_legacy_context.png"
        self.assertGreater(figure.stat().st_size, 100_000)


if __name__ == "__main__":
    unittest.main()
