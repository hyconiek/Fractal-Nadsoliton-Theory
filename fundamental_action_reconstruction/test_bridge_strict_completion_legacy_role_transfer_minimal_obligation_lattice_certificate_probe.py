import json
import subprocess
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "fundamental_action_reconstruction" / "scratch" / "bridge_strict_completion_legacy_role_transfer_minimal_obligation_lattice_certificate_probe.py"
REPORT_JSON = ROOT / "fundamental_action_reconstruction" / "scratch" / "bridge_strict_completion_legacy_role_transfer_minimal_obligation_lattice_certificate_report.json"
REPORT_MD = ROOT / "fundamental_action_reconstruction" / "scratch" / "bridge_strict_completion_legacy_role_transfer_minimal_obligation_lattice_certificate_report.md"


class LegacyRoleTransferMinimalObligationLatticeCertificateProbeTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        subprocess.run([sys.executable, str(SCRIPT)], check=True, cwd=ROOT)
        cls.payload = json.loads(REPORT_JSON.read_text(encoding="utf-8"))

    def test_result_kind_sources_and_atoms(self):
        payload = self.payload
        self.assertEqual(
            payload["result_kind"],
            "SCRATCH_STRICT_COMPLETION_LEGACY_ROLE_TRANSFER_MINIMAL_OBLIGATION_LATTICE_CERTIFICATE__NO_TRANSFER_THEOREM",
        )
        self.assertIn("minimal-role-transfer-obligation-lattice", payload["status"])
        self.assertEqual(set(payload["source_reports"]), {
            "role_transfer_pre_audit",
            "symbolic_cancellation",
            "component_gap_matrix",
            "legacy_bridge_guardrail",
        })
        self.assertEqual(payload["theorem_atoms"], [
            "alpha_geo_electroweak_role_theorem",
            "beta_tors_strict_role_theorem",
            "beta_power_hierarchy_successor_theorem",
            "chi11_selector_source_theorem",
        ])
        self.assertIn("minimal theorem", payload["grep_disambiguation"]["searched_terms"])
        self.assertIn("minimal missing theorem atoms", payload["grep_disambiguation"]["finding"])

    def test_lattice_rows_summary_and_cross_checks(self):
        payload = self.payload
        rows = payload["role_requirement_lattice_rows"]
        self.assertEqual(len(rows), 4)
        by_role = {row["role_id"]: row for row in rows}
        self.assertEqual(by_role["legacy_weak_mixing_angle"]["required_theorem_atoms"], ["alpha_geo_electroweak_role_theorem"])
        self.assertEqual(set(by_role["legacy_inverse_alpha_em"]["required_theorem_atoms"]), {"alpha_geo_electroweak_role_theorem", "beta_tors_strict_role_theorem"})
        self.assertEqual(set(by_role["legacy_beta_power_gravity_hierarchy"]["required_theorem_atoms"]), {"beta_power_hierarchy_successor_theorem", "beta_tors_strict_role_theorem"})
        self.assertEqual(set(by_role["legacy_torsion_to_chi11_orientation"]["required_theorem_atoms"]), {"beta_tors_strict_role_theorem", "chi11_selector_source_theorem"})
        self.assertTrue(all(not row["transfer_allowed_now"] for row in rows))

        coverage = {row["theorem_atom"]: row for row in payload["atom_coverage_rows"]}
        self.assertEqual(coverage["beta_tors_strict_role_theorem"]["covered_role_count"], 3)
        self.assertTrue(all(not row["exported_now"] for row in coverage.values()))

        summary = payload["role_transfer_minimal_obligation_summary"]
        self.assertEqual(summary["role_count"], 4)
        self.assertEqual(summary["theorem_atom_count"], 4)
        self.assertEqual(summary["total_subset_count_checked"], 16)
        self.assertTrue(summary["all_pre_audit_roles_loaded"])
        self.assertTrue(summary["all_roles_blocked_in_pre_audit"])
        self.assertTrue(summary["all_atoms_missing_now"])
        self.assertEqual(summary["global_minimal_obligation_count"], 1)
        self.assertEqual(summary["global_minimal_obligation_size"], 4)
        self.assertEqual(summary["global_minimal_obligation_sets"], [payload["theorem_atoms"]])
        self.assertTrue(summary["beta_tors_atom_is_shared_by_three_roles"])
        self.assertFalse(summary["role_transfer_theorem_exported"])
        self.assertFalse(summary["q2191_discharged"])
        self.assertFalse(summary["toe_closure_claimed"])
        self.assertTrue(all(payload["cross_checks"].values()))

    def test_proof_and_hard_limits(self):
        payload = self.payload
        proof = payload["proof_certificate"]
        self.assertIn("rg was used", proof["grep_step"])
        self.assertIn("All 16 subsets", proof["enumeration_step"])
        self.assertIn("zero atoms exported", proof["role_step"])
        self.assertIn("unique minimal obligation set", proof["global_step"])
        self.assertIn("beta_tors strict-role theorem is shared", proof["shared_beta_step"])
        self.assertIn("exports no role-transfer theorem", proof["scope_step"])

        hard_limits = "\n".join(payload["hard_limits"])
        self.assertIn("No theorem atom", hard_limits)
        self.assertIn("No legacy physical-role claim", hard_limits)
        self.assertIn("No alpha_geo, beta_tors, beta^N, or chi11", hard_limits)
        self.assertIn("No QW-2191 selector discharge", hard_limits)
        self.assertIn("No ToE closure", hard_limits)
        self.assertTrue(REPORT_MD.exists())


if __name__ == "__main__":
    unittest.main()
