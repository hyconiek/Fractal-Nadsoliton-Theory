from __future__ import annotations

import json
import subprocess
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "fundamental_action_reconstruction" / "scratch" / "bridge_completed_kernel_blocker_dependency_lattice_probe.py"
REPORT_JSON = ROOT / "fundamental_action_reconstruction" / "scratch" / "bridge_completed_kernel_blocker_dependency_lattice_report.json"
REPORT_MD = ROOT / "fundamental_action_reconstruction" / "scratch" / "bridge_completed_kernel_blocker_dependency_lattice_report.md"


class CompletedKernelBlockerDependencyLatticeProbeTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        subprocess.run([sys.executable, str(SCRIPT)], check=True, cwd=ROOT)
        cls.payload = json.loads(REPORT_JSON.read_text(encoding="utf-8"))

    def test_result_kind_lattice_and_cross_checks(self):
        payload = self.payload
        self.assertEqual(
            payload["result_kind"],
            "SCRATCH_COMPLETED_KERNEL_BLOCKER_DEPENDENCY_LATTICE__NO_FALSE_PASS",
        )
        self.assertEqual(payload["status"], "minimal-premise-lattice-computed-for-completed-strict-kernel-bridge-frontier")
        lattice = payload["finite_lattice"]
        self.assertEqual(lattice["premise_count"], 13)
        self.assertEqual(lattice["premise_subset_count"], 8192)
        self.assertEqual(lattice["exported_certificate_count"], 6)
        self.assertEqual(lattice["open_blocker_count"], 7)
        self.assertIn("rg was used", payload["repo_grep_nonduplication_note"])

        checks = payload["cross_report_checks"]
        self.assertEqual(checks["input_report_count"], 9)
        self.assertTrue(checks["completed_factorization_residual_pass"])
        self.assertTrue(checks["completed_factorization_current_live_kernel_is_strict"])
        self.assertTrue(checks["diagrams_reactivation_all_evidence_found"])
        self.assertEqual(checks["one_bit_frontier_name"], "T_beta_tors_orientation_exports_chi11_or_nonbridge")
        self.assertEqual(checks["torsion_opinion_overall_verdict"], "PARTLY_USEFUL_HEURISTIC_BUT_NOT_A_CURRENT_REPO_THEOREM")
        self.assertTrue(checks["reynolds_annihilator_zero"])
        self.assertEqual(checks["reynolds_chi11_rank"], 13)
        self.assertTrue(checks["local_puiseux_match_pass"])
        self.assertTrue(checks["eta_plus_one_match_pass"])
        self.assertLess(checks["measure_transport_balance_residual"], 1e-12)
        self.assertLess(checks["monotone_flow_output_matching_residual"], 1e-6)
        self.assertFalse(checks["cutoff_formal_inverse_global_admissible"])
        self.assertTrue(checks["global_admissibility_obstructed"])

    def test_premise_status_and_frontier_summary(self):
        payload = self.payload
        self.assertEqual(
            payload["exported_premises"],
            [
                "diagrams_legacy_carrier",
                "exact_completion_factorization",
                "local_puiseux_match",
                "eta_plus_one_puiseux_match",
                "measure_transport_identity",
                "monotone_flow_output_matching",
            ],
        )
        self.assertEqual(
            payload["open_blockers"],
            [
                "strict_transport_derivation",
                "global_z12_map_derivation",
                "orientation_chi11_source",
                "chi11_uniqueness",
                "reynolds_obstruction_escape",
                "role_transfer_theorem",
                "qw2191_selector_discharge",
            ],
        )
        status_by_premise = {row["premise"]: row["status"] for row in payload["premise_status_rows"]}
        self.assertEqual(status_by_premise["exact_completion_factorization"], "EXPORTED_CERTIFICATE")
        self.assertEqual(status_by_premise["orientation_chi11_source"], "OPEN_BLOCKER")
        self.assertEqual(status_by_premise["qw2191_selector_discharge"], "OPEN_BLOCKER")

        summary = payload["current_frontier_summary"]
        self.assertIn("K_strict_gate is the live full kernel", summary["current_live_kernel"])
        self.assertIn("orientation_chi11_source", summary["main_one_bit_blocker"])
        self.assertIn("strict_transport_derivation", summary["transport_blocker"])
        self.assertIn("role_transfer_theorem", summary["role_and_selector_blocker"])
        self.assertEqual(
            summary["remaining_for_theorem_level_bridge"],
            [
                "strict_transport_derivation",
                "global_z12_map_derivation",
                "orientation_chi11_source",
                "chi11_uniqueness",
                "reynolds_obstruction_escape",
                "role_transfer_theorem",
            ],
        )

    def test_minimal_premise_antichains(self):
        antichains = self.payload["minimal_premise_antichains"]
        self.assertEqual(antichains["finite_exact_completion_certificate"], [["exact_completion_factorization"]])
        self.assertEqual(antichains["guarded_completed_kernel_candidate"], [["diagrams_legacy_carrier", "exact_completion_factorization"]])
        self.assertEqual(antichains["local_two_channel_completion_candidate"], [["local_puiseux_match"]])
        self.assertEqual(antichains["eta_plus_one_local_completion_candidate"], [["local_puiseux_match", "eta_plus_one_puiseux_match"]])
        self.assertEqual(antichains["measure_transport_bookkeeping_candidate"], [["measure_transport_identity"]])
        self.assertEqual(antichains["monotone_output_matching_flow_candidate"], [["monotone_flow_output_matching"]])
        self.assertEqual(
            antichains["theorem_level_completed_legacy_to_strict_bridge"],
            [[
                "diagrams_legacy_carrier",
                "exact_completion_factorization",
                "strict_transport_derivation",
                "global_z12_map_derivation",
                "orientation_chi11_source",
                "chi11_uniqueness",
                "reynolds_obstruction_escape",
                "role_transfer_theorem",
            ]],
        )
        self.assertEqual(
            antichains["selector_closed_completed_kernel_toe_step"],
            [[
                "diagrams_legacy_carrier",
                "exact_completion_factorization",
                "strict_transport_derivation",
                "global_z12_map_derivation",
                "orientation_chi11_source",
                "chi11_uniqueness",
                "reynolds_obstruction_escape",
                "role_transfer_theorem",
                "qw2191_selector_discharge",
            ]],
        )
        self.assertEqual(antichains["strict_full_aut_internal_chi11_source"], [])

    def test_outcomes_proof_guardrails_and_markdown(self):
        rows = {row["outcome"]: row for row in self.payload["outcome_rows"]}
        self.assertTrue(rows["finite_exact_completion_certificate"]["currently_realized_by_loaded_reports"])
        self.assertTrue(rows["guarded_completed_kernel_candidate"]["currently_realized_by_loaded_reports"])
        self.assertFalse(rows["theorem_level_completed_legacy_to_strict_bridge"]["currently_realized_by_loaded_reports"])
        self.assertEqual(rows["theorem_level_completed_legacy_to_strict_bridge"]["minimal_premise_count"], 8)
        self.assertEqual(rows["selector_closed_completed_kernel_toe_step"]["minimal_premise_count"], 9)
        self.assertIsNone(rows["strict_full_aut_internal_chi11_source"]["minimal_premise_count"])

        proof = self.payload["exact_proof_certificate"]
        self.assertIn("2^13=8192", proof["finite_domain"])
        self.assertIn("already certificates", proof["separation_result"])
        self.assertIn("theorem-level completed legacy-to-strict bridge requires", proof["bridge_antichain"])
        self.assertIn("no minimal premise set", proof["full_aut_limit"])

        ontology = self.payload["ontology_guardrail"]
        self.assertIn("nadsoliton itself", ontology["allowed_reading"])
        self.assertIn("No separate informational layer", ontology["forbidden_reading"])

        hard_limits = "\n".join(self.payload["hard_limits"])
        self.assertIn("K_strict_gate is the current full form", hard_limits)
        self.assertIn("No unqualified identity K_legacy_ont == K_strict_gate", hard_limits)
        self.assertIn("No beta_tors -> chi_11 theorem", hard_limits)
        self.assertIn("No strict derivation", hard_limits)
        self.assertIn("No QW-2191 discharge", hard_limits)
        self.assertIn("No ToE closure", hard_limits)
        self.assertTrue(REPORT_MD.exists())


if __name__ == "__main__":
    unittest.main()
