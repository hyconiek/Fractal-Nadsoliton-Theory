from __future__ import annotations

import json
import subprocess
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "fundamental_action_reconstruction" / "scratch" / "bridge_strict_completion_s1_selector_margin_monotonicity_obstruction_certificate_probe.py"
REPORT_JSON = ROOT / "fundamental_action_reconstruction" / "scratch" / "bridge_strict_completion_s1_selector_margin_monotonicity_obstruction_certificate_report.json"
REPORT_MD = ROOT / "fundamental_action_reconstruction" / "scratch" / "bridge_strict_completion_s1_selector_margin_monotonicity_obstruction_certificate_report.md"


class StrictCompletionS1SelectorMarginMonotonicityObstructionCertificateProbeTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        subprocess.run([sys.executable, str(SCRIPT)], check=True, cwd=ROOT)
        cls.payload = json.loads(REPORT_JSON.read_text(encoding="utf-8"))

    def test_result_kind_sources_and_definition(self):
        payload = self.payload
        self.assertEqual(
            payload["result_kind"],
            "SCRATCH_STRICT_COMPLETION_S1_SELECTOR_MARGIN_MONOTONICITY_OBSTRUCTION_CERTIFICATE__FINITE_PREFLIGHT",
        )
        self.assertIn("s1-selector-margin", payload["status"])
        self.assertEqual(set(payload["source_reports"]), {
            "closure_plan_dependency_certificate",
            "p2275_corner_replay_passrate_floor",
            "p2277_adaptive_confidence_margin",
            "p2278_confidence_curve_sweep",
            "p2279_locked_confidence_profile_replay",
            "P1445_next_honest_step_packet",
        })
        self.assertIn("strict-local selector margin monotonicity witness", payload["grep_disambiguation"]["searched_terms"])
        definition = payload["s1_definition"]
        self.assertEqual(definition["source"], "P1445")
        self.assertEqual(definition["provider_class_audited"], "strict_nu_branch_group_policy_surrogate_confidence_route")
        self.assertIn("local surrogate preflight only", definition["scope"])
        self.assertEqual(len(definition["required_local_pass_conditions"]), 3)

    def test_finite_margin_table(self):
        payload = self.payload
        table = payload["finite_margin_table"]
        self.assertEqual(len(table["sweep_groups_by_delta"]), 4)
        self.assertEqual(len(table["floor_rows_by_risk"]), 4)
        self.assertEqual(table["locked_global_summary"]["worst_margin_to_target"], -0.99)
        for group in table["sweep_groups_by_delta"]:
            self.assertEqual(group["trial_multipliers"], [1, 2, 4])
            self.assertTrue(group["confidence_margin_decreases"])
            self.assertFalse(group["strictly_improves_worst_margin"])
            self.assertFalse(group["final_worst_margin_positive"])
        for row in table["floor_rows_by_risk"]:
            self.assertFalse(row["locked_meets_target"])
            self.assertLess(row["locked_margin_to_target"], 0)

    def test_summary_proof_blockers_and_hard_limits(self):
        payload = self.payload
        summary = payload["s1_obstruction_summary"]
        self.assertTrue(summary["closure_plan_recommends_s1"])
        self.assertTrue(summary["p1445_contains_s1_target"])
        self.assertTrue(summary["confidence_margins_decrease_with_budget"])
        self.assertFalse(summary["worst_selector_margins_strictly_improve_with_budget"])
        self.assertFalse(summary["worst_selector_margins_reach_positive"])
        self.assertFalse(summary["locked_replay_meets_all_targets"])
        self.assertEqual(summary["locked_worst_margin_to_target"], -0.99)
        self.assertFalse(summary["s1_witness_exported"])
        self.assertTrue(summary["s1_obstruction_certified"])

        blocker = "\n".join(payload["blocker_context"]["still_open"])
        self.assertIn("new provider class", blocker)
        self.assertIn("QW-2191", blocker)

        proof = payload["proof_certificate"]
        self.assertIn("Repo grep", proof["grep_step"])
        self.assertIn("never becomes positive", proof["sweep_step"])
        self.assertIn("worst_margin_to_target=-0.99", proof["locked_replay_step"])
        self.assertIn("no S1 witness is exported", proof["verdict_step"])
        self.assertIn("changing provider class", proof["next_step"])

        hard_limits = "\n".join(payload["hard_limits"])
        self.assertIn("No strict-local selector margin monotonicity witness", hard_limits)
        self.assertIn("No QW-2191 selector discharge", hard_limits)
        self.assertIn("No strict F_nadsoliton => L_SM + L_GR bridge", hard_limits)
        self.assertIn("No ToE closure", hard_limits)
        self.assertIn("No unqualified identity K_legacy_ont == K_strict_gate", hard_limits)
        self.assertTrue(REPORT_MD.exists())


if __name__ == "__main__":
    unittest.main()
