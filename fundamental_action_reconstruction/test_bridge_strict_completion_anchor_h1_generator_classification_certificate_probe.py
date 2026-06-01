from __future__ import annotations

import json
import subprocess
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "fundamental_action_reconstruction" / "scratch" / "bridge_strict_completion_anchor_h1_generator_classification_certificate_probe.py"
REPORT_JSON = ROOT / "fundamental_action_reconstruction" / "scratch" / "bridge_strict_completion_anchor_h1_generator_classification_certificate_report.json"
REPORT_MD = ROOT / "fundamental_action_reconstruction" / "scratch" / "bridge_strict_completion_anchor_h1_generator_classification_certificate_report.md"


class StrictCompletionAnchorH1GeneratorClassificationCertificateProbeTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        subprocess.run([sys.executable, str(SCRIPT)], check=True, cwd=ROOT)
        cls.payload = json.loads(REPORT_JSON.read_text(encoding="utf-8"))

    def test_result_kind_sources_and_type_ledger(self):
        payload = self.payload
        self.assertEqual(
            payload["result_kind"],
            "SCRATCH_STRICT_COMPLETION_ANCHOR_H1_GENERATOR_CLASSIFICATION_CERTIFICATE__GF2_TYPE_AUDIT",
        )
        self.assertIn("anchor-is-c0-gauge-fix", payload["status"])
        self.assertEqual(set(payload["source_reports"]), {
            "phase_sign_path_cohomology_triviality_certificate",
            "phase_sign_cycle_closure_obstruction_certificate",
            "phase_sign_z2_coboundary_certificate",
            "phase_sign_gf2_linear_system_certificate",
        })
        self.assertIn("cohomology generator", payload["grep_disambiguation"]["searched_terms"])

        ledger = payload["cohomology_type_ledger"]
        self.assertEqual(ledger["field"], "GF(2)")
        self.assertEqual(ledger["path_dim_C0"], 12)
        self.assertEqual(ledger["path_dim_C1"], 11)
        self.assertEqual(ledger["path_h1_dimension"], 0)
        self.assertIn("C0", ledger["path_anchor_object"])
        self.assertEqual(ledger["cycle_dim_C0"], 12)
        self.assertEqual(ledger["cycle_dim_C1"], 12)
        self.assertEqual(ledger["cycle_rank_delta"], 11)
        self.assertEqual(ledger["cycle_h1_dimension_dim_C1_minus_rank_delta"], 1)
        self.assertIn("C1 edge-cochain", ledger["cycle_generator_object"])
        self.assertEqual(ledger["image_size"], ledger["expected_image_size_2_rank"])

    def test_representatives_and_opinion_audit(self):
        payload = self.payload
        rows = {row["name"]: row for row in payload["representative_rows"]}
        audited = rows["audited_closed_phase_cochain"]
        self.assertEqual(audited["parity"], 0)
        self.assertTrue(audited["in_image_delta"])
        self.assertEqual(audited["h1_class"], "zero")

        closing = rows["closing_edge_generator_candidate"]
        self.assertEqual(closing["support"], ["11->0"])
        self.assertEqual(closing["parity"], 1)
        self.assertFalse(closing["in_image_delta"])
        self.assertEqual(closing["h1_class"], "nonzero_generator")

        first = rows["first_edge_generator_candidate"]
        self.assertEqual(first["parity"], 1)
        self.assertTrue(first["same_h1_class_as_closing_edge_generator"])

        even = rows["even_two_edge_boundary_candidate"]
        self.assertEqual(even["parity"], 0)
        self.assertTrue(even["in_image_delta"])
        self.assertEqual(even["h1_class"], "zero")

        verdicts = {row["claim"]: row["verdict"] for row in payload["opinion_audit_rows"]}
        self.assertIn("accepted_for_closed_cycle_only", verdicts.values())
        self.assertIn("rejected_type_mismatch", verdicts.values())
        self.assertIn("rejected_for_audited_zero_closure", verdicts.values())
        self.assertIn("accepted", verdicts.values())

    def test_summary_proof_and_hard_limits(self):
        payload = self.payload
        summary = payload["classification_summary"]
        self.assertTrue(summary["path_h1_zero"])
        self.assertTrue(summary["cycle_h1_one"])
        self.assertTrue(summary["cycle_image_equals_even_parity_kernel"])
        self.assertTrue(summary["closing_edge_odd_parity_generator"])
        self.assertTrue(summary["alternate_odd_edge_same_generator_class"])
        self.assertTrue(summary["audited_closed_cochain_is_exact_zero_h1_class"])
        self.assertTrue(summary["left_anchor_is_c0_gauge_fix_not_c1_generator"])
        self.assertTrue(summary["ai_opinion_as_stated_rejected"])
        self.assertTrue(summary["ai_opinion_weak_form_accepted"])
        self.assertTrue(summary["selector_source_remains_open"])
        self.assertTrue(summary["gf2_path_solution_unique_inherited"])

        proof = payload["proof_certificate"]
        self.assertIn("H^1(path;GF(2))=0", proof["path_step"])
        self.assertIn("C1/im(delta)", proof["cycle_step"])
        self.assertIn("odd parity represents the nonzero H^1 class", proof["parity_step"])
        self.assertIn("11->0", proof["generator_step"])
        self.assertIn("not a C1 edge cochain", proof["anchor_step"])
        self.assertIn("zero H^1 class", proof["audited_pattern_step"])
        self.assertIn("does not derive the selector/source bit", proof["theoretical_limit"])

        hard_limits = "\n".join(payload["hard_limits"])
        self.assertIn("No beta_tors -> chi_11 theorem", hard_limits)
        self.assertIn("No strict-core selector source", hard_limits)
        self.assertIn("No QW-2191 selector discharge", hard_limits)
        self.assertIn("No ToE closure", hard_limits)
        self.assertTrue(REPORT_MD.exists())


if __name__ == "__main__":
    unittest.main()
