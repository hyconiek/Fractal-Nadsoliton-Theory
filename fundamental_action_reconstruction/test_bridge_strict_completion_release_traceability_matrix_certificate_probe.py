import json
import subprocess
import sys
import unittest
from pathlib import Path


class StrictCompletionReleaseTraceabilityMatrixCertificateProbeTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.repo = Path(__file__).resolve().parents[1]
        cls.script = cls.repo / "fundamental_action_reconstruction" / "scratch" / "bridge_strict_completion_release_traceability_matrix_certificate_probe.py"
        subprocess.run([sys.executable, str(cls.script)], check=True, cwd=cls.repo)
        cls.report = cls.repo / "fundamental_action_reconstruction" / "scratch" / "bridge_strict_completion_release_traceability_matrix_certificate_report.json"
        cls.markdown = cls.repo / "fundamental_action_reconstruction" / "scratch" / "bridge_strict_completion_release_traceability_matrix_certificate_report.md"
        cls.payload = json.loads(cls.report.read_text(encoding="utf-8"))

    def test_result_kind_matrix_and_sources(self):
        payload = self.payload
        self.assertEqual(
            payload["result_kind"],
            "SCRATCH_STRICT_COMPLETION_RELEASE_TRACEABILITY_MATRIX_CERTIFICATE__GF2_SOURCE_TO_SECTION_AUDIT_NO_CLOSURE",
        )
        self.assertTrue(self.report.exists())
        self.assertTrue(self.markdown.exists())
        self.assertEqual(payload["matrix_columns"], ["D", "L", "R", "S"])
        self.assertEqual(len(payload["traceability_rows"]), 6)
        self.assertEqual(payload["column_coverage"], [4, 1, 3, 6])
        self.assertEqual(payload["row_coverage"], [2, 2, 2, 2, 3, 3])
        self.assertEqual(payload["gf2_rank"], 4)
        self.assertEqual(
            set(payload["source_reports"]),
            {"release_scaffold", "release_source_coverage", "chain_integrity", "frontier_low_weight", "anchor_h1"},
        )

    def test_rows_cross_checks_and_summary(self):
        payload = self.payload
        self.assertTrue(payload["all_cross_checks_pass"])
        for row in payload["traceability_rows"]:
            self.assertGreaterEqual(row["row_coverage"], 2)
            self.assertTrue(row["all_required_trace_snippets_present"])
            self.assertTrue(row["blocker_recorded"])

        checks = payload["cross_checks"]
        self.assertTrue(checks["traceability_doc_present"])
        self.assertTrue(checks["target_docs_present"])
        self.assertTrue(checks["required_trace_doc_snippets_present"])
        self.assertTrue(checks["matrix_shape_pass"])
        self.assertTrue(checks["column_coverage_pass"])
        self.assertTrue(checks["row_coverage_pass"])
        self.assertTrue(checks["gf2_rank_full_column_pass"])
        self.assertTrue(checks["release_source_coverage_inherited"])
        self.assertTrue(checks["release_scaffold_inherited"])
        self.assertTrue(checks["chain_integrity_inherited"])
        self.assertTrue(checks["selector_frontier_still_open"])
        self.assertTrue(checks["hard_limits_preserved"])

        summary = payload["traceability_summary"]
        self.assertTrue(summary["all_target_scaffold_columns_covered"])
        self.assertTrue(summary["all_source_rows_used_at_least_twice"])
        self.assertTrue(summary["target_columns_independent_over_gf2"])
        self.assertTrue(summary["source_coverage_and_scaffold_certificates_inherited"])
        self.assertTrue(summary["selector_frontier_still_open"])
        self.assertTrue(summary["no_bridge_theorem_claimed"])
        self.assertTrue(summary["no_full_eom_closure_claimed"])
        self.assertTrue(summary["no_role_transfer_claimed"])
        self.assertTrue(summary["no_qw2191_discharge"])
        self.assertTrue(summary["no_toe_closure"])

    def test_proof_and_limits(self):
        payload = self.payload
        proof = payload["proof_certificate"]
        self.assertIn("shape 6x4", proof["matrix_step"])
        self.assertIn("GF(2) rank 4", proof["matrix_step"])
        self.assertIn("legacy history", proof["source_step"])
        self.assertIn("strict Lagrangian/EOM", proof["source_step"])
        self.assertIn("chi11_selector_source", proof["frontier_step"])
        self.assertIn("no singleton/pair closes bridge", proof["frontier_step"])
        self.assertIn("traceability certificate only", proof["limit_step"])
        hard_limits = "\n".join(payload["hard_limits"])
        self.assertIn("No traceability edge is promoted to a theorem", hard_limits)
        self.assertIn("No bridge theorem", hard_limits)
        self.assertIn("No full tensor-resolved EOM closure", hard_limits)
        self.assertIn("No legacy physical-role transfer", hard_limits)
        self.assertIn("No QW-2191 selector discharge", hard_limits)
        self.assertIn("No ToE closure", hard_limits)


if __name__ == "__main__":
    unittest.main()
