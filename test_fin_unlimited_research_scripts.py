#!/usr/bin/env python3
"""Fast integrity tests for the unbounded P475/P485/P486/P487 package."""

from __future__ import annotations

import ast
from fractions import Fraction
import json
from pathlib import Path
import unittest

import fin_phase_exact_algebra as algebra


ROOT = Path(__file__).resolve().parent


class UnlimitedResearchScriptTests(unittest.TestCase):
    def test_p486_exact_premises_are_paid(self) -> None:
        result = json.loads(
            (ROOT / "FIN_Program_486_Results.json").read_text(encoding="utf-8")
        )
        self.assertTrue(result["orientation_premises_paid"])
        self.assertGreater(Fraction(result["det_C"]["lower"]), 0)
        self.assertLess(Fraction(result["active_pivot_determinant"]["upper"]), 0)
        self.assertLess(Fraction(result["rho_reference_coefficient"]["upper"]), 0)
        self.assertTrue(result["standard_sector_factorization_exact"])

    def test_exact_context_dimensions_and_selected_pivot(self) -> None:
        context = algebra.exact_context()
        tangent = algebra.tangent_consistency(context)
        pivot = algebra.best_active_pivot(context)
        self.assertEqual(len(context["active_equations"]), 12)
        self.assertEqual(len(context["tangent_equations"]), 6)
        self.assertEqual(len(tangent["consistency_polynomials"]), 5)
        self.assertEqual(tangent["reference"], 1)
        self.assertEqual(tuple(pivot["rows"]), (3, 5, 6))

    def test_unbounded_workers_have_no_timeout_or_rlimit_call(self) -> None:
        for name in (
            "fin_program_485_unlimited.py",
            "fin_program_487_unlimited.py",
            "fin_program_475_unlimited.py",
            "run_fin_p475_p485_p486_p487_unlimited.py",
        ):
            tree = ast.parse((ROOT / name).read_text(encoding="utf-8"))
            calls = [node for node in ast.walk(tree) if isinstance(node, ast.Call)]
            for call in calls:
                function = call.func
                attribute = function.attr if isinstance(function, ast.Attribute) else ""
                identifier = function.id if isinstance(function, ast.Name) else ""
                self.assertNotIn(attribute, {"setrlimit", "alarm"}, name)
                self.assertNotIn(identifier, {"timeout"}, name)
                self.assertFalse(any(keyword.arg == "timeout" for keyword in call.keywords), name)

    def test_pipeline_order_is_priority_and_memory_safe(self) -> None:
        source = (ROOT / "run_fin_p475_p485_p486_p487_unlimited.py").read_text(
            encoding="utf-8"
        )
        tree = ast.parse(source)
        namespace: dict[str, object] = {}
        stage_assignment = next(
            node for node in tree.body
            if isinstance(node, ast.Assign)
            and any(isinstance(target, ast.Name) and target.id == "STAGES"
                    for target in node.targets)
        )
        stages = ast.literal_eval(stage_assignment.value)
        self.assertEqual([row[0] for row in stages], [
            "P486", "P485", "P487", "P475-unlimited"
        ])
        self.assertIn("subprocess.run", source)
        self.assertNotIn("subprocess.Popen", source)


if __name__ == "__main__":
    unittest.main()
