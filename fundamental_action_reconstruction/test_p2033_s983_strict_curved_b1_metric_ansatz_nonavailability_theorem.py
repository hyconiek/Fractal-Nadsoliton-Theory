from __future__ import annotations

import json
import subprocess
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent
P1848 = ROOT / "p1848_s798_strict_gravity_componentwise_variation_and_counterterm_witness_checkpoint.py"
P2030 = ROOT / "p2030_s980_strict_tensor_resolved_projection_source_audit.py"
P2031 = ROOT / "p2031_s981_strict_b1_tensor_component_profile_table_scaffold.py"
P2032 = ROOT / "p2032_s982_strict_b1_metric_gauge_component_projection_rule_audit.py"
P2033 = ROOT / "p2033_s983_strict_curved_b1_metric_ansatz_nonavailability_theorem.py"
OUT2033 = ROOT / "generated" / "p2033_s983_strict_curved_b1_metric_ansatz_nonavailability_theorem.json"


class TestStrictCurvedB1MetricAnsatzNonavailabilityTheorem(unittest.TestCase):
    def test_p2033_proves_current_export_nonavailability_without_exporting_ansatz(self) -> None:
        subprocess.run([sys.executable, str(P1848)], check=True)
        subprocess.run([sys.executable, str(P2030)], check=True)
        subprocess.run([sys.executable, str(P2031)], check=True)
        subprocess.run([sys.executable, str(P2032)], check=True)
        subprocess.run([sys.executable, str(P2033)], check=True)
        data = json.loads(OUT2033.read_text(encoding="utf-8"))

        self.assertEqual(data["schema_version"], "p2033_s983_v1")
        self.assertEqual(data["status"], "OPEN_OBSTRUCTION_WITH_TRACE")
        self.assertEqual(
            data["result_kind"],
            "FORMAL_NONAVAILABILITY_OF_CURVED_B1_METRIC_ANSATZ_CURRENT_STRICT_EXPORTS",
        )
        self.assertEqual(data["branch_decision"]["decision"], "PROVE_NONAVAILABILITY_DO_NOT_EXPORT_MINIMAL_ANSATZ")
        self.assertFalse(data["branch_decision"]["minimal_curved_b1_ansatz_exported"])
        self.assertFalse(data["branch_decision"]["component_projection_rule_exported"])

        required = data["required_object"]
        self.assertEqual(
            required["object_id"],
            "minimal_curved_B1_metric_ansatz_and_component_projection_rule_v1",
        )
        self.assertEqual(required["current_export_state"], "NOT_EXPORTED")
        self.assertEqual(len(required["required_contract_fields"]), 6)
        self.assertEqual(len(required["missing_contract_fields"]), 6)
        self.assertIn("B1_background_metric_ansatz_g_munu_of_d", required["missing_contract_fields"])
        self.assertIn(
            "component_projection_operator_from_covariant_H_munu_templates",
            required["missing_contract_fields"],
        )
        self.assertIn("same_basis_one_loop_divergence_tensor_target", required["missing_contract_fields"])

        rows = data["near_miss_source_audit"]
        sources = {row["source"] for row in rows}
        self.assertTrue({"P1850", "P1870", "P1955", "P1958", "P1984", "P2032"}.issubset(sources))
        self.assertTrue(all(row["can_export_minimal_curved_B1_ansatz_now"] is False for row in rows))

        p1850 = next(row for row in rows if row["source"] == "P1850")
        self.assertIn("coefficient-only", p1850["exact_gap"])
        p1955 = next(row for row in rows if row["source"] == "P1955")
        self.assertIn("flat perturbative gauge-sector", p1955["exact_gap"])
        p1984 = next(row for row in rows if row["source"] == "P1984")
        self.assertIn("not H_00(d)", p1984["exact_gap"])

        theorem = data["formal_nonavailability_theorem"]
        self.assertTrue(theorem["export_state_nonavailability"])
        self.assertTrue(theorem["not_a_no_go_theorem"])
        self.assertIn("not available", theorem["statement"])
        self.assertIn("new premise/export", theorem["proof_trace"][-1])

        checks = data["gatekeeper_checks"]
        self.assertTrue(checks["p2031_table_remains_unfilled"])
        self.assertTrue(checks["p2032_all_required_fields_missing"])
        self.assertTrue(checks["all_near_miss_sources_available_or_audited"])
        self.assertTrue(checks["no_near_miss_source_can_export_ansatz"])
        self.assertTrue(checks["minimal_curved_b1_ansatz_not_exported"])
        self.assertTrue(checks["component_projection_rule_not_exported"])
        self.assertTrue(checks["nonavailability_theorem_passed"])
        self.assertTrue(checks["no_tensor_component_profile_filled"])
        self.assertTrue(checks["no_tensor_projection_claimed"])
        self.assertTrue(checks["no_independent_a_GB_claimed"])
        self.assertTrue(checks["no_toe_closure_claimed"])

        not_licensed = data["theorem_scope"]["not_licensed"]
        self.assertIn("choosing a new g_munu(d) ansatz without marking it as a new premise/export", not_licensed)
        self.assertIn("using B1 coefficient formulas as metric components", not_licensed)
        self.assertIn("using ADM lapse cancellation as H_00(d)", not_licensed)


if __name__ == "__main__":
    unittest.main()
