from __future__ import annotations

import json
import subprocess
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent
P1848 = ROOT / "p1848_s798_strict_gravity_componentwise_variation_and_counterterm_witness_checkpoint.py"
P2027 = ROOT / "p2027_s977_strict_b1_rank3_gb_null_adaptive_quadrature_witness.py"
P2028 = ROOT / "p2028_s978_strict_b1_gb_quotient_counterterm_identifiability_theorem.py"
P2029 = ROOT / "p2029_s979_strict_task1_renormalization_quotient_ledger_update.py"
P2030 = ROOT / "p2030_s980_strict_tensor_resolved_projection_source_audit.py"
P2031 = ROOT / "p2031_s981_strict_b1_tensor_component_profile_table_scaffold.py"
P2032 = ROOT / "p2032_s982_strict_b1_metric_gauge_component_projection_rule_audit.py"
P2033 = ROOT / "p2033_s983_strict_curved_b1_metric_ansatz_nonavailability_theorem.py"
P2034 = ROOT / "p2034_s984_strict_task1_quotient_only_renormalization_theorem.py"
OUT2034 = ROOT / "generated" / "p2034_s984_strict_task1_quotient_only_renormalization_theorem.json"


class TestStrictTask1QuotientOnlyRenormalizationTheorem(unittest.TestCase):
    def test_p2034_licenses_only_local_scalar_b1_quotient_renormalization(self) -> None:
        for script in (P1848, P2027, P2028, P2029, P2030, P2031, P2032, P2033, P2034):
            subprocess.run([sys.executable, str(script)], check=True)

        data = json.loads(OUT2034.read_text(encoding="utf-8"))

        self.assertEqual(data["schema_version"], "p2034_s984_v1")
        self.assertEqual(data["status"], "OPEN_OBSTRUCTION_WITH_TRACE")
        self.assertEqual(
            data["result_kind"],
            "PASS_LOCAL_B1_QUOTIENT_ONLY_RENORMALIZATION_THEOREM__GLOBAL_TENSOR_OPEN",
        )
        self.assertEqual(data["local_verdict"], "PASS_QUOTIENT_ONLY_RENORMALIZATION_WITH_TRACE")

        self.assertFalse(data["grep_audit_decision"]["p2034_already_done"])
        self.assertEqual(
            data["professor_decision"]["decision"],
            "PIVOT_TO_QUOTIENT_ONLY_RENORMALIZATION_DO_NOT_ADD_CURVED_B1_PREMISE",
        )
        self.assertEqual(data["professor_decision"]["rejected_branch_for_now"], "new-premise curved B1 metric ansatz contract")

        theorem = data["quotient_renormalization_theorem"]
        self.assertEqual(theorem["status"], "PASS_WITH_TRACE")
        self.assertEqual(theorem["quotient_space"], "R^4 / span(1,-4,1,-1)")
        self.assertEqual(theorem["null_direction_R2_Ric2_Riem2_GB"], [1.0, -4.0, 1.0, -1.0])
        self.assertIn("T(a_R2,a_Ric2,a_Riem2,a_GB)", theorem["quotient_map"])
        self.assertIn("not physical measurements", theorem["section_warning"])

        patch = data["task1_p2029_patch_effect"]
        self.assertTrue(patch["p2029_explicit_quotient_option_present"])
        self.assertTrue(patch["local_quotient_option_now_licensed"])
        self.assertTrue(patch["does_not_satisfy_tensor_or_global_items"])
        self.assertIn("background-global transport", patch["remaining_open_items"][3])

        numeric = data["numeric_trace"]
        self.assertTrue(numeric["scalar_b1_numerics_stable"])
        self.assertLessEqual(numeric["tau_family_max_full_system_residual_linf"], numeric["numerical_tolerance"])
        self.assertLessEqual(numeric["tau_family_max_quotient_gap_linf"], numeric["numerical_tolerance"])

        checks = data["gatekeeper_checks"]
        self.assertTrue(checks["rank3_quotient_input_pass"])
        self.assertTrue(checks["scalar_b1_numerics_stable"])
        self.assertTrue(checks["tensor_branch_blocked_current_exports"])
        self.assertTrue(checks["p2031_table_unfilled"])
        self.assertTrue(checks["no_new_curved_b1_metric_premise_added"])
        self.assertTrue(checks["renormalization_licensed_only_in_quotient"])
        self.assertTrue(checks["independent_a_GB_not_identified"])
        self.assertTrue(checks["tensor_projection_not_claimed"])
        self.assertTrue(checks["no_tensor_component_profile_filled"])
        self.assertTrue(checks["no_background_globalization_claimed"])
        self.assertTrue(checks["no_toe_closure_claimed"])

        not_licensed = data["theorem_scope"]["not_licensed"]
        self.assertIn("new curved B1 metric ansatz", not_licensed)
        self.assertIn("any tensor component profile H_00/H_11/H_22/H_33", not_licensed)
        self.assertIn("four-channel counterterm uniqueness", not_licensed)
        self.assertIn("ToE closure", not_licensed)


if __name__ == "__main__":
    unittest.main()
