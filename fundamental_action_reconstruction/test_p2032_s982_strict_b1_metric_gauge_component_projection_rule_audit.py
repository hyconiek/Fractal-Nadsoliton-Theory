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
OUT2032 = ROOT / "generated" / "p2032_s982_strict_b1_metric_gauge_component_projection_rule_audit.json"


class TestStrictB1MetricGaugeComponentProjectionRuleAudit(unittest.TestCase):
    def test_p2032_blocks_component_projection_rule_without_curved_b1_ansatz(self) -> None:
        subprocess.run([sys.executable, str(P1848)], check=True)
        subprocess.run([sys.executable, str(P2030)], check=True)
        subprocess.run([sys.executable, str(P2031)], check=True)
        subprocess.run([sys.executable, str(P2032)], check=True)
        data = json.loads(OUT2032.read_text(encoding="utf-8"))

        self.assertEqual(data["schema_version"], "p2032_s982_v1")
        self.assertEqual(data["status"], "OPEN_OBSTRUCTION_WITH_TRACE")
        self.assertEqual(data["result_kind"], "OPEN_B1_METRIC_GAUGE_COMPONENT_PROJECTION_RULE_MISSING_WITH_TRACE")
        self.assertEqual(data["audited_rule"]["rule_id"], "strict_B1_metric_gauge_component_projection_rule_v1")
        self.assertFalse(data["audited_rule"]["rule_ready"])

        required = data["audited_rule"]["required_rule_fields"]
        missing = data["audited_rule"]["missing_rule_fields"]
        self.assertEqual(len(required), 6)
        self.assertEqual(len(missing), 6)
        self.assertIn("B1_background_metric_ansatz_g_munu_of_d", missing)
        self.assertIn("component_projection_operator_from_covariant_H_munu_templates", missing)
        self.assertIn("same_basis_one_loop_divergence_tensor_target", missing)

        rows = data["candidate_source_rows"]
        self.assertEqual(len(rows), 6)
        self.assertTrue(all(row["can_supply_rule"] is False for row in rows))

        positive = data["positive_but_insufficient_evidence"]
        self.assertTrue(positive["p1848_covariant_templates_present"])
        self.assertTrue(positive["p1950_scalar_b1_projection_present"])
        self.assertTrue(positive["p1954_metric_perturbation_missing_recorded"])
        self.assertTrue(positive["p1958_flat_tangent_metric_present"])
        self.assertTrue(positive["p1984_gb_lapse_cancellation_present"])
        self.assertTrue(positive["p2031_table_scaffold_missing_all_entries"])

        checks = data["gatekeeper_checks"]
        self.assertTrue(checks["candidate_sources_audited"])
        self.assertTrue(checks["no_candidate_source_can_supply_rule"])
        self.assertTrue(checks["required_rule_fields_all_missing"])
        self.assertTrue(checks["flat_tangent_metric_not_promoted"])
        self.assertTrue(checks["gb_lapse_not_promoted_to_H00"])
        self.assertTrue(checks["p2031_remains_unfilled"])
        self.assertFalse(checks["rule_ready"])
        self.assertTrue(checks["no_tensor_projection_claimed"])
        self.assertTrue(checks["no_toe_closure_claimed"])

        self.assertIn("deriving H_00/H_ii from flat tangent-patch eta alone", data["theorem_scope"]["not_licensed"])
        self.assertIn("deriving H_00 from ADM lapse Euler cancellation", data["theorem_scope"]["not_licensed"])


if __name__ == "__main__":
    unittest.main()
