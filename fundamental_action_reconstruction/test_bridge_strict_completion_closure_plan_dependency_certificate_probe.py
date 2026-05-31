from __future__ import annotations

import json
import subprocess
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "fundamental_action_reconstruction" / "scratch" / "bridge_strict_completion_closure_plan_dependency_certificate_probe.py"
REPORT_JSON = ROOT / "fundamental_action_reconstruction" / "scratch" / "bridge_strict_completion_closure_plan_dependency_certificate_report.json"
REPORT_MD = ROOT / "fundamental_action_reconstruction" / "scratch" / "bridge_strict_completion_closure_plan_dependency_certificate_report.md"


class StrictCompletionClosurePlanDependencyCertificateProbeTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        subprocess.run([sys.executable, str(SCRIPT)], check=True, cwd=ROOT)
        cls.payload = json.loads(REPORT_JSON.read_text(encoding="utf-8"))

    def test_result_kind_sources_and_grep(self):
        payload = self.payload
        self.assertEqual(
            payload["result_kind"],
            "SCRATCH_STRICT_COMPLETION_CLOSURE_PLAN_DEPENDENCY_CERTIFICATE__REPO_FRONTIER_MATRIX_NO_CLOSURE_CLAIM",
        )
        self.assertIn("closure-plan-dependency", payload["status"])
        self.assertIn("strict_completion_chain_integrity_report", payload["source_reports"])
        self.assertIn("P1445_next_honest_step", payload["source_reports"])
        self.assertIn("N679_selector_boundary", payload["source_reports"])
        self.assertIn("domknięcie mostu", payload["grep_disambiguation"]["searched_terms"])
        self.assertIn("no single strict-completion closure-plan dependency matrix", payload["grep_disambiguation"]["finding"])

    def test_dependency_matrix_and_plan(self):
        payload = self.payload
        certificate = payload["dependency_certificate"]
        self.assertEqual(certificate["recommended_next_step"], "legacy_to_strict_completion_bridge_guardrail")
        self.assertEqual(certificate["recommended_next_step_source"], "S2_strategic_reorientation")
        self.assertEqual(certificate["plan_steps_in_topological_order"][0], "finite_strict_completion_ledger")
        self.assertEqual(certificate["plan_steps_in_topological_order"][-1], "strict_core_ToE_closure")
        self.assertEqual(len(certificate["dependency_matrix_rows"]), 7)
        self.assertEqual(certificate["next_open_steps_in_order"][0], "legacy_to_strict_completion_bridge_guardrail")
        self.assertEqual(certificate["next_selector_step_after_bridge_guardrail"], "strict_local_selector_margin_monotonicity_witness")
        for row in certificate["dependency_matrix_rows"]:
            self.assertEqual(len(row["dependency_bits_in_plan_order"]), 7)

        plan = payload["closure_plan_for_bridge_and_theory"]
        self.assertEqual([row["order"] for row in plan], [1, 2, 3, 4, 5, 6])
        self.assertIn("legacy -> strict", plan[0]["step"])
        self.assertIn("selector-margin witness", plan[1]["step"])
        self.assertIn("QW-2191", plan[3]["step"])
        self.assertIn("F_nadsoliton => L_SM + L_GR", plan[4]["step"])
        self.assertIn("ToE closure", plan[5]["step"])

    def test_summary_proof_and_hard_limits(self):
        payload = self.payload
        checks = payload["source_keyword_checks"]
        self.assertTrue(all(checks.values()))
        summary = payload["closure_plan_summary"]
        self.assertTrue(summary["chain_ledger_currently_cross_consistent"])
        self.assertTrue(summary["all_source_keyword_checks_pass"])
        self.assertTrue(summary["dependency_matrix_is_triangular_in_plan_order"])
        self.assertTrue(summary["recommended_next_step_is_legacy_strict_bridge_guardrail"])
        self.assertTrue(summary["S1_selector_margin_remains_next_selector_subproblem"])
        self.assertTrue(summary["role_transfer_audit_required"])
        self.assertTrue(summary["qw2191_or_orientation_remains_hard_blocker"])
        self.assertTrue(summary["toe_closure_not_claimed"])
        self.assertTrue(summary["legacy_strict_identity_not_used"])

        proof = payload["proof_certificate"]
        self.assertIn("Repo grep", proof["grep_step"])
        self.assertIn("finite strict-completion ledger", proof["ledger_step"])
        self.assertIn("triangular", proof["matrix_step"])
        self.assertIn("legacy->strict completion-bridge guardrail", proof["next_step"])
        self.assertIn("S1 remains", proof["next_step"])
        self.assertIn("QW-2191", proof["blocker_step"])
        self.assertIn("does not complete the legacy->strict bridge", proof["theoretical_limit"])
        self.assertIn("prove S1", proof["theoretical_limit"])

        hard_limits = "\n".join(payload["hard_limits"])
        self.assertIn("No unqualified identity K_legacy_ont == K_strict_gate", hard_limits)
        self.assertIn("No strict-local selector margin witness", hard_limits)
        self.assertIn("No QW-2191 selector discharge", hard_limits)
        self.assertIn("No strict F_nadsoliton => L_SM + L_GR bridge", hard_limits)
        self.assertIn("No ToE closure", hard_limits)
        self.assertTrue(REPORT_MD.exists())


if __name__ == "__main__":
    unittest.main()
