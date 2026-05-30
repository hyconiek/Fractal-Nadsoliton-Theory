from __future__ import annotations

import json
import subprocess
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "fundamental_action_reconstruction" / "scratch" / "bridge_legacy_kernel_reactivation_from_diagrams_candidate_probe.py"
REPORT_JSON = ROOT / "fundamental_action_reconstruction" / "scratch" / "bridge_legacy_kernel_reactivation_from_diagrams_candidate_report.json"
REPORT_MD = ROOT / "fundamental_action_reconstruction" / "scratch" / "bridge_legacy_kernel_reactivation_from_diagrams_candidate_report.md"


class LegacyKernelReactivationFromDiagramsCandidateProbeTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        subprocess.run([sys.executable, str(SCRIPT)], check=True, cwd=ROOT)
        cls.payload = json.loads(REPORT_JSON.read_text(encoding="utf-8"))

    def test_result_kind_status_and_cross_checks(self):
        payload = self.payload
        self.assertEqual(
            payload["result_kind"],
            "SCRATCH_LEGACY_KERNEL_REACTIVATION_FROM_DIAGRAMS_CANDIDATE__THEOREM_TARGET_NOT_CLOSURE",
        )
        self.assertEqual(
            payload["status"],
            "legacy-kernel-restored-as-guarded-bridge-theorem-target-not-silent-identity",
        )
        self.assertIn("2026-05-30", payload["author_session_update"])
        checks = payload["source_cross_checks"]
        self.assertTrue(checks["diagrams_evidence_all_found"])
        self.assertTrue(checks["k1_kernel_split_note_present"])
        self.assertTrue(checks["s2_retirement_recorded_but_reinterpreted_by_current_author_instruction"])
        self.assertTrue(checks["f253_future_bridge_target_exists"])
        self.assertTrue(checks["chi11_reduced_bit_supported"])
        self.assertTrue(checks["chi11_meta_uniqueness_still_open"])
        self.assertTrue(checks["reynolds_annihilator_still_blocks_full_aut_source"])
        self.assertEqual(checks["conditional_balanced_ledger"], [2, 2, 2, 1, 1])
        self.assertEqual(checks["target_q_power"], "256/243")
        self.assertEqual(checks["target_eta"], "9/5")

    def test_diagrams_evidence_and_reactivation_decision(self):
        evidence = {row["evidence_id"]: row for row in self.payload["diagrams_evidence_rows"]}
        self.assertEqual(
            set(evidence),
            {
                "k_total_four_mechanisms",
                "effective_kernel_compression",
                "torsion_cosine_fingerprint",
                "hyperbolic_damping_fingerprint",
                "beta_tors_sensitivity",
            },
        )
        self.assertTrue(all(row["found"] for row in evidence.values()))
        self.assertIn("nadsoliton-characteristic carrier", evidence["effective_kernel_compression"]["bridge_relevance"])
        self.assertIn("chi_11 source", evidence["torsion_cosine_fingerprint"]["bridge_relevance"])

        decision = self.payload["reactivation_decision"]
        self.assertIn("positive bridge theorem-target", decision["current_author_update_reading"])
        self.assertIn("candidate completed/renormalized", decision["new_working_position"])
        self.assertIn("Do not state K_legacy_ont == K_strict_gate", decision["not_allowed"])

    def test_bridge_slots_one_bit_frontier_and_requirements(self):
        slots = {row["slot"]: row for row in self.payload["bridge_slot_table"]}
        self.assertEqual(
            set(slots),
            {
                "legacy_characteristic_carrier",
                "amplitude_alpha_geo",
                "torsion_orientation_bit",
                "character_uniqueness",
                "reynolds_obstruction_escape",
                "eta_compression_pipeline",
            },
        )
        self.assertEqual(slots["legacy_characteristic_carrier"]["current_status"], "REACTIVATED_AS_THEOREM_TARGET_BY_AUTHOR_INSTRUCTION")
        self.assertEqual(slots["torsion_orientation_bit"]["current_status"], "MOST_IMPORTANT_ONE_BIT_BRIDGE_TARGET")
        self.assertEqual(slots["torsion_orientation_bit"]["missing_object"], "T_beta_tors_orientation_exports_chi11_or_nonbridge")
        self.assertEqual(slots["character_uniqueness"]["current_status"], "OPEN_UNIQUENESS_BLOCKER")
        self.assertEqual(slots["reynolds_obstruction_escape"]["current_status"], "MUST_BE_NON_FULL_AUT_OR_EXPLICITLY_ORIENTED")
        self.assertEqual(slots["eta_compression_pipeline"]["current_status"], "CONDITIONAL_DOWNSTREAM_CERTIFICATES_EXIST")

        frontier = self.payload["one_bit_frontier"]
        self.assertEqual(frontier["frontier_name"], "T_beta_tors_orientation_exports_chi11_or_nonbridge")
        self.assertIn("non-full-Aut orientation source", frontier["positive_branch"])
        self.assertIn("nonbridge theorem", frontier["negative_branch"])
        self.assertIn("DIAGRAMS supplies", frontier["why_this_is_now_the_main_bridge_subproblem"])

        requirements = "\n".join(self.payload["theorem_target_requirements"])
        self.assertIn("legacy torsion/orientation observable", requirements)
        self.assertIn("uniquely chi_11", requirements)
        self.assertIn("beta_tors*d to beta*d^eta", requirements)
        self.assertIn("completed K_legacy_ont", requirements)

    def test_proof_ontology_guardrails_and_markdown(self):
        proof = self.payload["exact_proof_certificate"]
        self.assertIn("restoring legacy K(d) as an active bridge theorem-target", proof["positive_upgrade"])
        self.assertIn("beta_tors/K_tors/topology orientation -> chi_11", proof["one_bit_focus"])
        self.assertIn("orientation map", proof["still_open"])
        self.assertIn("candidate completed/renormalized legacy kernel", proof["guarded_language"])

        ontology = self.payload["ontology_guardrail"]
        self.assertIn("nadsoliton itself", ontology["allowed_reading"])
        self.assertIn("No separate informational layer", ontology["forbidden_reading"])

        hard_limits = "\n".join(self.payload["hard_limits"])
        self.assertIn("No unqualified identity K_legacy_ont == K_strict_gate", hard_limits)
        self.assertIn("No beta_tors -> chi_11 theorem", hard_limits)
        self.assertIn("No legacy physical-role transfer", hard_limits)
        self.assertIn("No theorem derives eta=9/5", hard_limits)
        self.assertIn("No QW-2191 discharge", hard_limits)
        self.assertIn("No ToE closure", hard_limits)
        self.assertTrue(REPORT_MD.exists())


if __name__ == "__main__":
    unittest.main()
