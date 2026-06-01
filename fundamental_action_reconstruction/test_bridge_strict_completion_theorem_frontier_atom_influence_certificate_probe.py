import json
import subprocess
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "fundamental_action_reconstruction" / "scratch" / "bridge_strict_completion_theorem_frontier_atom_influence_certificate_probe.py"
REPORT_JSON = ROOT / "fundamental_action_reconstruction" / "scratch" / "bridge_strict_completion_theorem_frontier_atom_influence_certificate_report.json"
REPORT_MD = ROOT / "fundamental_action_reconstruction" / "scratch" / "bridge_strict_completion_theorem_frontier_atom_influence_certificate_report.md"


class TheoremFrontierAtomInfluenceCertificateProbeTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        subprocess.run([sys.executable, str(SCRIPT)], check=True, cwd=ROOT)
        cls.payload = json.loads(REPORT_JSON.read_text(encoding="utf-8"))

    def test_result_kind_sources_and_scope(self):
        payload = self.payload
        self.assertEqual(
            payload["result_kind"],
            "SCRATCH_STRICT_COMPLETION_THEOREM_FRONTIER_ATOM_INFLUENCE_CERTIFICATE__NO_THEOREM_EXPORT",
        )
        self.assertIn("boolean-swing-criticality", payload["status"])
        self.assertEqual(set(payload["source_reports"]), {
            "theorem_frontier_truth_table",
            "theorem_frontier_cut",
            "role_transfer_minimal_obligation_lattice",
        })
        self.assertEqual(len(payload["open_atoms"]), 7)
        self.assertEqual(payload["targets"], [
            "bridge_theorem_level_closure",
            "role_transfer_theorem_level_closure",
            "selector_qw2191_closure",
            "toe_closure",
        ])
        self.assertIn("atom influence", payload["grep_disambiguation"]["searched_terms"])
        self.assertIn("finite influence audit", payload["grep_disambiguation"]["finding"])

    def test_influence_rows_summary_and_cross_checks(self):
        payload = self.payload
        rows = {row["atom"]: row for row in payload["atom_influence_rows"]}
        self.assertEqual(len(rows), 7)
        self.assertEqual(payload["atom_influence_order"][0], "chi11_selector_source")
        self.assertEqual(rows["chi11_selector_source"]["selector_qw2191_critical_count"], 64)
        self.assertEqual(rows["chi11_selector_source"]["role_transfer_critical_count"], 8)
        self.assertEqual(rows["chi11_selector_source"]["toe_critical_count"], 1)
        self.assertEqual(rows["chi11_selector_source"]["total_critical_count"], 73)
        for atom in [
            "strict_dynamical_source_for_A_P_D",
            "strict_phase_frequency_source",
            "strict_damping_beta_eta_source",
        ]:
            self.assertEqual(rows[atom]["bridge_critical_count"], 16)
            self.assertEqual(rows[atom]["toe_critical_count"], 1)
            self.assertEqual(rows[atom]["total_critical_count"], 17)
        for atom in [
            "alpha_geo_electroweak_role_theorem",
            "beta_tors_strict_role_theorem",
            "beta_power_hierarchy_successor_theorem",
        ]:
            self.assertEqual(rows[atom]["role_transfer_critical_count"], 8)
            self.assertEqual(rows[atom]["toe_critical_count"], 1)
            self.assertEqual(rows[atom]["total_critical_count"], 9)

        summary = payload["theorem_frontier_atom_influence_summary"]
        self.assertEqual(summary["open_atom_count"], 7)
        self.assertEqual(summary["target_count"], 4)
        self.assertEqual(summary["critical_pair_universe_per_atom_target"], 64)
        self.assertEqual(summary["total_critical_count_all_atoms"], 151)
        self.assertEqual(summary["top_influence_atoms"], ["chi11_selector_source"])
        self.assertEqual(summary["top_influence_total_critical_count"], 73)
        self.assertTrue(summary["chi11_selector_source_is_unique_top_influence"])
        self.assertTrue(summary["bridge_source_atoms_tie_at_17"])
        self.assertTrue(summary["role_only_atoms_tie_at_9"])
        self.assertTrue(summary["each_atom_is_toe_critical_once"])
        self.assertTrue(summary["truth_table_current_assignment_closes_no_target"])
        self.assertTrue(summary["role_lattice_atoms_remain_missing"])
        self.assertFalse(summary["bridge_theorem_exported"])
        self.assertFalse(summary["role_transfer_theorem_exported"])
        self.assertFalse(summary["selector_closure_exported"])
        self.assertFalse(summary["toe_closure_claimed"])
        self.assertTrue(all(payload["cross_checks"].values()))

    def test_proof_and_hard_limits(self):
        payload = self.payload
        proof = payload["proof_certificate"]
        self.assertIn("rg was used", proof["grep_step"])
        self.assertIn("flipping a to true", proof["definition_step"])
        self.assertIn("7*4*64", proof["enumeration_step"])
        self.assertIn("total critical count 73", proof["influence_step"])
        self.assertIn("total critical count 17", proof["bridge_step"])
        self.assertIn("total critical count 9", proof["role_step"])
        self.assertIn("exports no atom", proof["scope_step"])

        hard_limits = "\n".join(payload["hard_limits"])
        self.assertIn("No influence score", hard_limits)
        self.assertIn("No missing theorem atom", hard_limits)
        self.assertIn("No bridge theorem", hard_limits)
        self.assertIn("No QW-2191 selector discharge", hard_limits)
        self.assertTrue(REPORT_MD.exists())


if __name__ == "__main__":
    unittest.main()
