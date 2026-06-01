import json
import subprocess
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "fundamental_action_reconstruction" / "scratch" / "bridge_strict_completion_legacy_role_transfer_pre_audit_certificate_probe.py"
REPORT_JSON = ROOT / "fundamental_action_reconstruction" / "scratch" / "bridge_strict_completion_legacy_role_transfer_pre_audit_certificate_report.json"
REPORT_MD = ROOT / "fundamental_action_reconstruction" / "scratch" / "bridge_strict_completion_legacy_role_transfer_pre_audit_certificate_report.md"


class LegacyRoleTransferPreAuditCertificateProbeTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        subprocess.run([sys.executable, str(SCRIPT)], check=True, cwd=ROOT)
        cls.payload = json.loads(REPORT_JSON.read_text(encoding="utf-8"))

    def test_result_kind_sources_and_columns(self):
        payload = self.payload
        self.assertEqual(
            payload["result_kind"],
            "SCRATCH_STRICT_COMPLETION_LEGACY_ROLE_TRANSFER_PRE_AUDIT_CERTIFICATE__NO_TRANSFER_THEOREM",
        )
        self.assertIn("legacy-physical-role-claims-audited", payload["status"])
        self.assertEqual(set(payload["source_reports"]), {
            "legacy_bridge_guardrail",
            "component_gap_matrix",
            "symbolic_cancellation",
            "amplitude_scalar_normalization",
            "damping_parameter_identifiability",
            "anchor_h1_classification",
        })
        self.assertEqual(set(payload["source_specs"]), {"s2_priority_packet", "t15_bridge_theorem_spec"})
        self.assertIn("role-transfer", payload["grep_disambiguation"]["searched_terms"])
        self.assertIn("blocks every listed legacy role", payload["grep_disambiguation"]["finding"])
        self.assertIn("depends_alpha_geo", payload["role_dependency_columns"])
        self.assertIn("depends_beta_tors", payload["role_dependency_columns"])

    def test_role_rows_matrix_summary_and_cross_checks(self):
        payload = self.payload
        rows = payload["role_transfer_rows"]
        self.assertEqual(len(rows), 4)
        self.assertEqual({row["role_id"] for row in rows}, {
            "legacy_weak_mixing_angle",
            "legacy_inverse_alpha_em",
            "legacy_beta_power_gravity_hierarchy",
            "legacy_torsion_to_chi11_orientation",
        })
        self.assertTrue(all(row["current_verdict"] == "blocked_not_transferred" for row in rows))
        self.assertTrue(all(not row["role_transfer_allowed_now"] for row in rows))
        self.assertEqual(len(payload["role_dependency_matrix_gf2"]), 4)
        self.assertTrue(any(row["depends_alpha_geo"] for row in rows))
        self.assertTrue(any(row["depends_beta_tors"] for row in rows))
        self.assertTrue(any(row["depends_chi11_or_orientation"] for row in rows))

        summary = payload["role_transfer_pre_audit_summary"]
        self.assertEqual(summary["role_claim_count"], 4)
        self.assertEqual(summary["dependency_matrix_rank_gf2"], 4)
        self.assertTrue(summary["all_roles_blocked_now"])
        self.assertEqual(summary["roles_transferred_now"], 0)
        self.assertTrue(summary["s2_lists_required_role_transfer_claims"])
        self.assertTrue(summary["t15_keeps_role_transfer_separate"])
        self.assertTrue(summary["guardrail_requires_role_transfer_audit"])
        self.assertTrue(summary["symbolic_bridge_is_formula_only"])
        self.assertTrue(summary["component_gap_blocks_role_transfer"])
        self.assertTrue(summary["alpha_role_source_missing"])
        self.assertTrue(summary["beta_tors_role_source_missing"])
        self.assertTrue(summary["chi11_source_missing"])
        self.assertFalse(summary["role_transfer_theorem_exported"])
        self.assertFalse(summary["q2191_discharged"])
        self.assertFalse(summary["toe_closure_claimed"])
        self.assertTrue(all(payload["cross_checks"].values()))

    def test_proof_and_hard_limits(self):
        payload = self.payload
        proof = payload["proof_certificate"]
        self.assertIn("rg was used", proof["grep_step"])
        self.assertIn("dependency matrix of rank 4", proof["matrix_step"])
        self.assertIn("sin^2(theta_W)=alpha_geo/12 is blocked", proof["alpha_step"])
        self.assertIn("alpha_EM^-1 and beta^N hierarchy claims are blocked", proof["beta_step"])
        self.assertIn("beta_tors -> chi_11 is blocked", proof["chi11_step"])
        self.assertIn("exports zero role-transfer permissions", proof["scope_step"])

        hard_limits = "\n".join(payload["hard_limits"])
        self.assertIn("No legacy physical-role claim", hard_limits)
        self.assertIn("No alpha_geo electroweak role theorem", hard_limits)
        self.assertIn("No beta_tors electromagnetic/gravity hierarchy", hard_limits)
        self.assertIn("No beta_tors -> chi_11 theorem", hard_limits)
        self.assertIn("No QW-2191 selector discharge", hard_limits)
        self.assertIn("No ToE closure", hard_limits)
        self.assertTrue(REPORT_MD.exists())


if __name__ == "__main__":
    unittest.main()
