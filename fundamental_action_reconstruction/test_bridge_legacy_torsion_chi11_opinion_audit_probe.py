from __future__ import annotations

import json
import subprocess
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "fundamental_action_reconstruction" / "scratch" / "bridge_legacy_torsion_chi11_opinion_audit_probe.py"
REPORT_JSON = ROOT / "fundamental_action_reconstruction" / "scratch" / "bridge_legacy_torsion_chi11_opinion_audit_report.json"
REPORT_MD = ROOT / "fundamental_action_reconstruction" / "scratch" / "bridge_legacy_torsion_chi11_opinion_audit_report.md"


class LegacyTorsionChi11OpinionAuditProbeTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        subprocess.run([sys.executable, str(SCRIPT)], check=True, cwd=ROOT)
        cls.payload = json.loads(REPORT_JSON.read_text(encoding="utf-8"))

    def test_result_kind_and_cross_checks(self):
        payload = self.payload
        self.assertEqual(
            payload["result_kind"],
            "SCRATCH_LEGACY_TORSION_CHI11_OPINION_AUDIT__CANDIDATE_NOT_THEOREM",
        )
        self.assertEqual(payload["status"], "ai-opinion-classified-as-interesting-but-overstated-bridge-hypothesis")
        checks = payload["finite_model_cross_checks"]
        self.assertEqual(checks["ring"], "Z_12")
        self.assertEqual(checks["support_count"], 792)
        self.assertEqual(checks["automorphism_units"], [1, 5, 7, 11])
        self.assertEqual(checks["chi11_character_count"], 4)
        self.assertEqual(checks["reynolds_rank_from_reynolds_probe"], 25)
        self.assertEqual(checks["chi11_rank_from_reynolds_probe"], 13)
        self.assertEqual(checks["strict_full_aut_internal_chi11_antichain"], [])
        self.assertEqual(checks["target_q_power"], "256/243")
        self.assertEqual(checks["target_eta"], "9/5")
        self.assertEqual(checks["balanced_ledger"], [2, 2, 2, 1, 1])

    def test_opinion_verdicts(self):
        audit = self.payload["opinion_audit"]
        self.assertEqual(audit["overall_verdict"], "PARTLY_USEFUL_HEURISTIC_BUT_NOT_A_CURRENT_REPO_THEOREM")
        claims = {row["claim_id"]: row for row in audit["claims"]}
        self.assertEqual(len(claims), 8)
        self.assertEqual(claims["legacy_beta_tors_exists_as_torsion_damping"]["verdict"], "SUPPORTED_AS_LEGACY_ARCHIVAL_FACT")
        self.assertTrue(claims["legacy_beta_tors_exists_as_torsion_damping"]["evidence"][0])
        self.assertEqual(claims["strict_kernel_has_different_beta_eta_form"]["verdict"], "SUPPORTED_STRUCTURAL_SPLIT_ONLY")
        self.assertEqual(claims["torsion_implies_parity_chirality_breaking"]["verdict"], "NOT_EXPORTED_BY_CURRENT_REPO")
        self.assertEqual(claims["chi11_is_missing_unit_character_bit"]["verdict"], "SUPPORTED_IN_REDUCED_ALGEBRAIC_AUDITS")
        self.assertEqual(claims["beta_tors_collapses_to_chi11"]["verdict"], "CANDIDATE_BRIDGE_HYPOTHESIS_NOT_THEOREM")
        self.assertEqual(claims["legacy_kernel_bridge_is_now_the_main_bridge"]["verdict"], "CONFLICTS_WITH_CURRENT_S2_STRATEGY_IF_ASSERTED_AS_FACT")
        self.assertEqual(claims["chi11_triggers_eta_and_balanced_ledger"]["verdict"], "OVERSTATED_CONDITIONAL_DOWNSTREAM_ONLY")
        self.assertEqual(claims["full_aut_invariant_source_exports_chi11"]["verdict"], "REFUTED_FOR_CURRENT_FULL_AUT_SUPPORT_DATA")
        self.assertTrue(claims["full_aut_invariant_source_exports_chi11"]["evidence"][0])
        self.assertTrue(claims["full_aut_invariant_source_exports_chi11"]["evidence"][2])

    def test_requirements_proof_guardrails_and_markdown(self):
        requirements = {row["requirement"]: row for row in self.payload["bridge_hypothesis_requirements"]}
        self.assertEqual(set(requirements), {"orientation_map", "role_transfer_control", "eta_pipeline_link", "full_aut_obstruction_escape"})
        self.assertIn("why the kernel is {1,11}", requirements["orientation_map"]["needed_theorem"])
        self.assertIn("retires", requirements["role_transfer_control"]["current_status"])
        self.assertIn("conditional", requirements["eta_pipeline_link"]["current_status"])
        self.assertIn("Reynolds annihilator", requirements["full_aut_obstruction_escape"]["needed_theorem"])

        recommended = self.payload["recommended_next_honest_step"]
        self.assertIn("non-strict theorem-target spec", recommended["if_pursuing_this_intuition"])
        self.assertIn("strict-only", recommended["if_following_current_S2_priority"])

        proof = self.payload["exact_proof_certificate"]
        self.assertIn("beta_tors as a legacy torsion/damping datum", proof["supported_part"])
        self.assertIn("No loaded report exports beta_tors -> chi_11", proof["unsupported_part"])
        self.assertIn("Full-Aut Reynolds averaging annihilates chi_11", proof["current_obstruction"])
        self.assertIn("retired to archival status", proof["s2_boundary"])

        ontology = self.payload["ontology_guardrail"]
        self.assertIn("nadsoliton itself", ontology["allowed_reading"])
        self.assertIn("No separate informational layer", ontology["forbidden_reading"])

        hard_limits = "\n".join(self.payload["hard_limits"])
        self.assertIn("No identity K_legacy_ont == K_strict_gate", hard_limits)
        self.assertIn("No beta_tors -> chi_11", hard_limits)
        self.assertIn("No legacy physical-role transfer", hard_limits)
        self.assertIn("No theorem derives eta=9/5", hard_limits)
        self.assertIn("No QW-2191 discharge", hard_limits)
        self.assertIn("No ToE closure", hard_limits)
        self.assertTrue(REPORT_MD.exists())


if __name__ == "__main__":
    unittest.main()
