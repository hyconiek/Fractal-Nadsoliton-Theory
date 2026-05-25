from __future__ import annotations

import json
import subprocess
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent
P2027 = ROOT / "p2027_s977_strict_b1_rank3_gb_null_adaptive_quadrature_witness.py"
P2028 = ROOT / "p2028_s978_strict_b1_gb_quotient_counterterm_identifiability_theorem.py"
P2029 = ROOT / "p2029_s979_strict_task1_renormalization_quotient_ledger_update.py"
P2030 = ROOT / "p2030_s980_strict_tensor_resolved_projection_source_audit.py"
OUT2029 = ROOT / "generated" / "p2029_s979_strict_task1_renormalization_quotient_ledger_update.json"
OUT2030 = ROOT / "generated" / "p2030_s980_strict_tensor_resolved_projection_source_audit.json"


class TestTask1TensorSourceDecisionChain(unittest.TestCase):
    def test_p2029_task1_ledger_preserves_quotient_not_four_channel_closure(self) -> None:
        subprocess.run([sys.executable, str(P2027)], check=True)
        subprocess.run([sys.executable, str(P2028)], check=True)
        subprocess.run([sys.executable, str(P2029)], check=True)
        data = json.loads(OUT2029.read_text(encoding="utf-8"))

        self.assertEqual(data["schema_version"], "p2029_s979_v1")
        self.assertEqual(data["status"], "OPEN_OBSTRUCTION_WITH_TRACE")
        self.assertEqual(data["result_kind"], "TASK1_LEDGER_UPDATED_TO_RANK3_QUOTIENT_PASS_TENSOR_EXTENSION_OPEN")
        self.assertTrue(data["task1_quotient_status"]["quotient_class_pass"])
        self.assertTrue(data["task1_quotient_status"]["adaptive_scalar_projection_pass"])
        self.assertFalse(data["task1_quotient_status"]["independent_a_GB_identified"])
        self.assertFalse(data["task1_quotient_status"]["tensor_resolved_projection_exported"])
        self.assertIn("tensor-resolved B1 curvature-operator projection", data["task1_ledger_patch"]["replacement_missing"][0])
        self.assertIn("minimum-norm four-coefficient vector", data["task1_ledger_patch"]["do_not_use_anymore"][1])
        self.assertEqual(data["local_readiness_score_0_1"], 0.6)
        self.assertTrue(data["gatekeeper_checks"]["no_toe_closure_claimed"])

    def test_p2030_tensor_source_audit_blocks_projection_without_component_table(self) -> None:
        subprocess.run([sys.executable, str(P2029)], check=True)
        subprocess.run([sys.executable, str(P2030)], check=True)
        data = json.loads(OUT2030.read_text(encoding="utf-8"))

        self.assertEqual(data["schema_version"], "p2030_s980_v1")
        self.assertEqual(data["status"], "OPEN_OBSTRUCTION_WITH_TRACE")
        self.assertEqual(data["result_kind"], "OPEN_TENSOR_RESOLVED_SOURCE_GAP_WITH_TRACE")
        self.assertFalse(data["tensor_projection_ready"])
        self.assertIn("tensor_component_profile_table", data["blocking_requirements"])
        self.assertIn("tensor_component_Gram_rule", data["blocking_requirements"])
        self.assertIn("same_basis_divergence_tensor_target", data["blocking_requirements"])
        self.assertTrue(data["known_positive_evidence"]["gb_lapse_cancels"])
        self.assertTrue(data["gatekeeper_checks"]["p1848_has_only_templates_not_components"])
        self.assertTrue(data["gatekeeper_checks"]["adm_is_lapse_not_full_tensor"])
        self.assertTrue(data["gatekeeper_checks"]["no_tensor_projection_claimed"])
        self.assertIn("ADM lapse cancellation", data["false_pass_guard"])
        self.assertIn("tensor-resolved Gram rank or coefficient identification", data["theorem_scope"]["not_licensed"])


if __name__ == "__main__":
    unittest.main()
