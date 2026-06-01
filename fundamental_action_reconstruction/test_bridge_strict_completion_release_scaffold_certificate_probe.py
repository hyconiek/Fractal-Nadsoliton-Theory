import json
import subprocess
import sys
import unittest
from pathlib import Path


class StrictCompletionReleaseScaffoldCertificateProbeTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.repo = Path(__file__).resolve().parents[1]
        cls.script = cls.repo / "fundamental_action_reconstruction" / "scratch" / "bridge_strict_completion_release_scaffold_certificate_probe.py"
        subprocess.run([sys.executable, str(cls.script)], check=True, cwd=cls.repo)
        cls.report = cls.repo / "fundamental_action_reconstruction" / "scratch" / "bridge_strict_completion_release_scaffold_certificate_report.json"
        cls.markdown = cls.repo / "fundamental_action_reconstruction" / "scratch" / "bridge_strict_completion_release_scaffold_certificate_report.md"
        cls.payload = json.loads(cls.report.read_text(encoding="utf-8"))

    def test_report_shape_and_sources(self):
        payload = self.payload
        self.assertEqual(
            payload["result_kind"],
            "SCRATCH_STRICT_COMPLETION_RELEASE_SCAFFOLD_CERTIFICATE__THREE_MD_FILES_AUDITED_NO_CLOSURE",
        )
        self.assertTrue(self.report.exists())
        self.assertTrue(self.markdown.exists())
        self.assertEqual(
            set(payload["scaffold_docs"]),
            {
                "strict_kernel_transformation_diagrams",
                "strict_kernel_lagrangian_eom_draft",
                "strict_kernel_role_transfer_audit_draft",
            },
        )
        self.assertEqual(
            set(payload["source_reports"]),
            {
                "chain_integrity",
                "finite_bridge_assembly",
                "symbolic_cancellation",
                "role_transfer_pre_audit",
                "role_transfer_minimal_obligation_lattice",
                "theorem_frontier_low_weight_extension",
                "anchor_h1_generator_classification",
                "strict_lagrangian_symbolic_export_p1866",
                "strict_schematic_eom_spectrum_p2315",
            },
        )

    def test_doc_checks_pass(self):
        payload = self.payload
        for doc_check in payload["doc_checks"].values():
            self.assertTrue(doc_check["exists"])
            self.assertTrue(doc_check["all_required_snippets_present"])
            self.assertEqual(
                doc_check["matched_required_snippet_count"],
                doc_check["required_snippet_count"],
            )
            self.assertGreater(doc_check["line_count"], 30)

    def test_cross_checks_and_summary(self):
        payload = self.payload
        self.assertTrue(payload["all_cross_checks_pass"])
        checks = payload["cross_checks"]
        self.assertTrue(checks["three_scaffold_files_present"])
        self.assertTrue(checks["strict_diagram_has_bridge_history"])
        self.assertTrue(checks["lagrangian_eom_draft_has_p1866_exports"])
        self.assertTrue(checks["role_transfer_audit_blocks_all_legacy_roles"])
        self.assertTrue(checks["chain_integrity_component_witnesses_present"])
        self.assertTrue(checks["finite_bridge_reconstructs_without_role_transfer"])
        self.assertTrue(checks["symbolic_cancellation_formula_only"])
        self.assertTrue(checks["role_transfer_pre_audit_blocks_all_roles"])
        self.assertTrue(checks["role_transfer_lattice_has_four_atom_obligation"])
        self.assertTrue(checks["low_weight_frontier_keeps_bridge_role_toe_open"])
        self.assertTrue(checks["anchor_h1_type_audit_preserves_selector_gap"])
        self.assertTrue(checks["p1866_symbolic_lagrangian_export_loaded"])
        self.assertTrue(checks["p2315_schematic_eom_keeps_qw2191_open"])
        self.assertTrue(checks["no_identity_bridge_claimed"])
        self.assertTrue(checks["no_role_transfer_theorem_claimed"])
        self.assertTrue(checks["no_qw2191_discharge"])
        self.assertTrue(checks["no_toe_closure"])

        summary = payload["release_scaffold_summary"]
        self.assertTrue(summary["strict_kernel_diagram_scaffold_ready"])
        self.assertTrue(summary["strict_lagrangian_eom_scaffold_ready"])
        self.assertTrue(summary["strict_role_transfer_audit_scaffold_ready"])
        self.assertTrue(summary["finite_bridge_ledger_inherited"])
        self.assertTrue(summary["symbolic_lagrangian_eom_exports_inherited"])
        self.assertTrue(summary["all_legacy_roles_blocked_pending_theorems"])
        self.assertTrue(summary["no_identity_bridge_claimed"])
        self.assertTrue(summary["no_role_transfer_theorem_claimed"])
        self.assertTrue(summary["no_qw2191_discharge"])
        self.assertTrue(summary["no_toe_closure"])

    def test_proof_certificate_and_hard_limits(self):
        payload = self.payload
        proof = payload["proof_certificate"]
        self.assertIn("DIAGRAMS_STRICT_KERNEL_TRANSFORMATION.md", proof["strict_diagram_step"])
        self.assertIn("K_legacy_ont", proof["strict_diagram_step"])
        self.assertIn("K_strict_gate", proof["strict_diagram_step"])
        self.assertIn("P1866", proof["lagrangian_eom_step"])
        self.assertIn("P2315", proof["lagrangian_eom_step"])
        self.assertIn("sin^2(theta_W)", proof["role_transfer_step"])
        self.assertIn("beta_tors->chi_11", proof["role_transfer_step"])
        self.assertIn("no legacy role is transferred", proof["role_transfer_step"])
        self.assertIn("release-build package", proof["ledger_step"])
        hard_limits = "\n".join(payload["hard_limits"])
        self.assertIn("No unqualified identity", hard_limits)
        self.assertIn("No legacy physical-role transfer theorem", hard_limits)
        self.assertIn("No QW-2191 selector discharge", hard_limits)
        self.assertIn("No ToE closure", hard_limits)


if __name__ == "__main__":
    unittest.main()
