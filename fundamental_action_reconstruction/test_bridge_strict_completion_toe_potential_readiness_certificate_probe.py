import json
import subprocess
import sys
import unittest
from pathlib import Path


class StrictCompletionToEPotentialReadinessCertificateProbeTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.repo = Path(__file__).resolve().parents[1]
        cls.script = cls.repo / "fundamental_action_reconstruction" / "scratch" / "bridge_strict_completion_toe_potential_readiness_certificate_probe.py"
        subprocess.run([sys.executable, str(cls.script)], check=True, cwd=cls.repo)
        cls.report = cls.repo / "fundamental_action_reconstruction" / "scratch" / "bridge_strict_completion_toe_potential_readiness_certificate_report.json"
        cls.markdown = cls.repo / "fundamental_action_reconstruction" / "scratch" / "bridge_strict_completion_toe_potential_readiness_certificate_report.md"
        cls.doc = cls.repo / "fundamental_action_reconstruction" / "STRICT_KERNEL_TOE_POTENTIAL_AUDIT.md"
        cls.payload = json.loads(cls.report.read_text(encoding="utf-8"))

    def test_result_kind_sources_and_frontier_rows(self):
        payload = self.payload
        self.assertEqual(
            payload["result_kind"],
            "SCRATCH_STRICT_COMPLETION_TOE_POTENTIAL_READINESS_CERTIFICATE__PROFESSORIAL_FINITE_AUDIT_NO_CLOSURE",
        )
        self.assertTrue(self.report.exists())
        self.assertTrue(self.markdown.exists())
        self.assertTrue(self.doc.exists())
        self.assertEqual(len(payload["open_atoms"]), 7)
        self.assertEqual(
            set(payload["source_reports"]),
            {
                "chain_integrity",
                "frontier_cut",
                "frontier_truth_table",
                "frontier_atom_influence",
                "frontier_target_signature_lattice",
                "frontier_low_weight_extension",
                "role_transfer_lattice",
                "release_source_coverage",
                "release_traceability_matrix",
            },
        )
        row_values = {row["metric"]: row["value"] for row in payload["frontier_quantitative_rows"]}
        self.assertEqual(row_values["open_atom_count"], 7)
        self.assertEqual(row_values["truth_assignment_count"], 128)
        self.assertEqual(row_values["toe_satisfying_assignment_count"], 1)
        self.assertEqual(row_values["toe_readiness_fraction"], "1/128")
        self.assertEqual(row_values["toe_minimal_set_size"], 7)
        self.assertEqual(row_values["reachable_target_signature_count"], 6)
        self.assertEqual(row_values["chi11_total_critical_count"], 73)

    def test_cross_checks_and_summary(self):
        payload = self.payload
        self.assertEqual(payload["status"], "PASS")
        self.assertTrue(payload["all_cross_checks_pass"])
        self.assertTrue(all(payload["cross_checks"].values()))
        checks = payload["cross_checks"]
        self.assertTrue(checks["toe_audit_doc_present"])
        self.assertTrue(checks["toe_audit_doc_required_snippets_present"])
        self.assertTrue(checks["source_reports_present"])
        self.assertTrue(checks["chain_integrity_inherited"])
        self.assertTrue(checks["open_atom_set_pass"])
        self.assertTrue(checks["truth_table_counts_pass"])
        self.assertTrue(checks["toe_minimal_weight_pass"])
        self.assertTrue(checks["target_signature_lattice_pass"])
        self.assertTrue(checks["low_weight_no_go_pass"])
        self.assertTrue(checks["chi11_top_bottleneck_pass"])
        self.assertTrue(checks["role_transfer_still_blocked"])
        self.assertTrue(checks["release_coverage_inherited"])
        self.assertTrue(checks["traceability_matrix_inherited"])
        self.assertTrue(checks["hard_limits_preserved"])

        summary = payload["toe_potential_readiness_summary"]
        self.assertTrue(summary["professorial_toe_potential_doc_certified"])
        self.assertTrue(summary["finite_frontier_board_certified"])
        self.assertTrue(summary["toe_requires_all_7_open_atoms"])
        self.assertEqual(summary["toe_readiness_fraction"], "1/128")
        self.assertEqual(summary["reachable_target_signatures"], ["0000", "0010", "0110", "1000", "1010", "1111"])
        self.assertTrue(summary["low_weight_extensions_do_not_close_toe"])
        self.assertTrue(summary["chi11_selector_source_is_top_bottleneck"])
        self.assertTrue(summary["role_transfer_remains_blocked"])
        self.assertTrue(summary["traceability_and_source_coverage_inherited"])
        self.assertFalse(summary["toe_closure_claimed"])

    def test_proof_and_limits(self):
        payload = self.payload
        proof = payload["proof_certificate"]
        self.assertIn("professor-level assessment", proof["professorial_assessment_step"])
        self.assertIn("2^7=128 assignments", proof["truth_table_step"])
        self.assertIn("1/128", proof["truth_table_step"])
        self.assertIn("all seven open atoms", proof["truth_table_step"])
        self.assertIn("6 reachable signatures out of 16", proof["signature_lattice_step"])
        self.assertIn("all 7 singleton and 21 pair extensions", proof["low_weight_step"])
        self.assertIn("unique top bottleneck", proof["bottleneck_step"])
        self.assertIn("total critical count 73", proof["bottleneck_step"])
        self.assertIn("finite readiness bookkeeping only", proof["limit_step"])
        hard_limits = "\n".join(payload["hard_limits"])
        self.assertIn("No identity bridge", hard_limits)
        self.assertIn("No strict dynamical source theorem", hard_limits)
        self.assertIn("No legacy physical-role transfer", hard_limits)
        self.assertIn("No beta_tors -> chi11 theorem", hard_limits)
        self.assertIn("No QW-2191 selector discharge", hard_limits)
        self.assertIn("No full tensor-resolved Lagrangian/EOM theorem", hard_limits)
        self.assertIn("No ToE closure", hard_limits)


if __name__ == "__main__":
    unittest.main()
