import json
import subprocess
import sys
import unittest
from pathlib import Path


class StrictCompletionReleaseSourceCoverageCertificateProbeTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.repo = Path(__file__).resolve().parents[1]
        cls.script = cls.repo / "fundamental_action_reconstruction" / "scratch" / "bridge_strict_completion_release_source_coverage_certificate_probe.py"
        subprocess.run([sys.executable, str(cls.script)], check=True, cwd=cls.repo)
        cls.report = cls.repo / "fundamental_action_reconstruction" / "scratch" / "bridge_strict_completion_release_source_coverage_certificate_report.json"
        cls.markdown = cls.repo / "fundamental_action_reconstruction" / "scratch" / "bridge_strict_completion_release_source_coverage_certificate_report.md"
        cls.payload = json.loads(cls.report.read_text(encoding="utf-8"))

    def test_result_kind_and_family_rows(self):
        payload = self.payload
        self.assertEqual(
            payload["result_kind"],
            "SCRATCH_STRICT_COMPLETION_RELEASE_SOURCE_COVERAGE_CERTIFICATE__REPO_GREP_PROVENANCE_NO_CLOSURE",
        )
        self.assertTrue(self.report.exists())
        self.assertTrue(self.markdown.exists())
        self.assertGreater(payload["scan_file_count"], 100)
        self.assertEqual(
            {row["family"] for row in payload["source_family_rows"]},
            {
                "legacy_kernel_history",
                "strict_lagrangian_eom",
                "finite_bridge_ledger",
                "role_transfer_selector_boundaries",
            },
        )
        for row in payload["source_family_rows"]:
            self.assertGreater(row["grep_like_hit_count"], 0)
            self.assertTrue(row["required_files_present"])
            self.assertTrue(row["required_files_mentioned_in_doc"])
            self.assertTrue(row["top_hits"])

    def test_cross_checks_and_summary(self):
        payload = self.payload
        self.assertTrue(payload["all_cross_checks_pass"])
        checks = payload["cross_checks"]
        self.assertTrue(checks["source_coverage_doc_present"])
        self.assertTrue(checks["all_required_doc_snippets_present"])
        self.assertTrue(checks["all_source_families_have_hits"])
        self.assertTrue(checks["all_required_family_files_present"])
        self.assertTrue(checks["all_required_family_files_mentioned"])
        self.assertTrue(checks["p2316_repo_grep_audit_available"])
        self.assertTrue(checks["p2316_keeps_full_task3_theorem_open"])
        self.assertTrue(checks["p1866_symbolic_export_still_obstruction"])
        self.assertTrue(checks["release_scaffold_already_passes"])
        self.assertTrue(checks["no_identity_role_selector_toe_closure"])

        summary = payload["coverage_summary"]
        self.assertTrue(summary["legacy_history_covered"])
        self.assertTrue(summary["strict_lagrangian_eom_sources_covered"])
        self.assertTrue(summary["finite_bridge_ledger_covered"])
        self.assertTrue(summary["role_transfer_selector_boundaries_covered"])
        self.assertTrue(summary["p2316_full_task3_theorem_still_open"])
        self.assertTrue(summary["release_scaffold_nonduplicating_source_map_ready"])
        self.assertTrue(summary["no_bridge_theorem_claimed"])
        self.assertTrue(summary["no_role_transfer_claimed"])
        self.assertTrue(summary["no_qw2191_discharge"])
        self.assertTrue(summary["no_toe_closure"])

    def test_proof_and_limits(self):
        payload = self.payload
        proof = payload["proof_certificate"]
        self.assertIn("broad rg pattern", proof["grep_step"])
        self.assertIn("four source families", proof["grep_step"])
        self.assertIn("DIAGRAMS_KERNEL_TRANSFORMATION.md", proof["coverage_step"])
        self.assertIn("P1622", proof["coverage_step"])
        self.assertIn("P1866", proof["coverage_step"])
        self.assertIn("P2315", proof["coverage_step"])
        self.assertIn("P2316", proof["coverage_step"])
        self.assertIn("N87/N103", proof["coverage_step"])
        self.assertIn("scaffolds", proof["nonduplication_step"])
        self.assertIn("OPEN_OBSTRUCTION_WITH_TRACE", proof["limit_step"])
        hard_limits = "\n".join(payload["hard_limits"])
        self.assertIn("No source-coverage hit is promoted to a theorem", hard_limits)
        self.assertIn("No bridge theorem", hard_limits)
        self.assertIn("No full tensor-resolved EOM closure", hard_limits)
        self.assertIn("No QW-2191 selector discharge", hard_limits)
        self.assertIn("No ToE closure", hard_limits)


if __name__ == "__main__":
    unittest.main()
