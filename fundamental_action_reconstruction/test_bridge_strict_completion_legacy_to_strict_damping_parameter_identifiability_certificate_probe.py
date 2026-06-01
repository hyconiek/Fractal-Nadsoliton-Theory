from __future__ import annotations

import json
import subprocess
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "fundamental_action_reconstruction" / "scratch" / "bridge_strict_completion_legacy_to_strict_damping_parameter_identifiability_certificate_probe.py"
REPORT_JSON = ROOT / "fundamental_action_reconstruction" / "scratch" / "bridge_strict_completion_legacy_to_strict_damping_parameter_identifiability_certificate_report.json"
REPORT_MD = ROOT / "fundamental_action_reconstruction" / "scratch" / "bridge_strict_completion_legacy_to_strict_damping_parameter_identifiability_certificate_report.md"


class LegacyToStrictDampingParameterIdentifiabilityCertificateProbeTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        subprocess.run([sys.executable, str(SCRIPT)], check=True, cwd=ROOT)
        cls.payload = json.loads(REPORT_JSON.read_text(encoding="utf-8"))

    def test_result_kind_sources_and_definition(self):
        payload = self.payload
        self.assertEqual(
            payload["result_kind"],
            "SCRATCH_STRICT_COMPLETION_LEGACY_TO_STRICT_DAMPING_PARAMETER_IDENTIFIABILITY_CERTIFICATE__NO_SOURCE_THEOREM",
        )
        self.assertIn("strict-beta-eta-finitely-identifiable", payload["status"])
        self.assertEqual(set(payload["source_reports"]), {
            "damping_compression_separation",
            "damping_exact_rational",
            "finite_diagonal_completion_map",
            "component_gap_matrix",
            "legacy_bridge_guardrail",
        })
        self.assertIn("eta identifiability", payload["grep_disambiguation"]["searched_terms"])
        self.assertIn("finite beta/eta identifiability", payload["grep_disambiguation"]["finding"])

        definition = payload["parameter_identification_definition"]
        self.assertIn("S(d)=1+beta*d^eta", definition["strict_denominator"])
        self.assertEqual(definition["beta_recovery"], "beta=S(1)-1=1")
        self.assertIn("log(S(d)-1)", definition["eta_recovery"])
        self.assertIn("(S(d)-1)^5=d^9", definition["finite_cover_identity"])

    def test_rows_summary_and_cross_checks(self):
        payload = self.payload
        rows = payload["parameter_recovery_rows"]
        self.assertEqual(len(rows), 11)
        self.assertEqual([row["d"] for row in rows], list(range(1, 12)))
        self.assertTrue(all(abs(row["beta_from_d_if_eta_fixed_9_over_5"] - 1.0) <= 1e-12 for row in rows))
        self.assertTrue(all(row["exact_cover_identity"] == "(S(d)-1)^5=d^9" for row in rows))
        eta_rows = [row for row in rows if row["eta_from_log_if_beta_fixed_1"] is not None]
        self.assertTrue(all(abs(row["eta_from_log_residual_vs_9_over_5"]) <= 1e-12 for row in eta_rows))
        gammas = [row["required_linear_gamma_for_this_node"] for row in rows]
        self.assertTrue(all(a < b for a, b in zip(gammas, gammas[1:])))

        summary = payload["damping_parameter_identifiability_summary"]
        self.assertEqual(summary["positive_domain"], list(range(1, 12)))
        self.assertTrue(summary["beta_identified_from_d1"])
        self.assertEqual(summary["identified_beta"], "1")
        self.assertEqual(summary["identified_eta"], "9/5")
        self.assertTrue(summary["all_beta_recoveries_equal_1_given_eta"])
        self.assertTrue(summary["all_eta_log_recoveries_equal_9_over_5_given_beta"])
        self.assertTrue(summary["exact_fifth_power_cover_identity_recorded"])
        self.assertEqual(summary["candidate_grid_p_le_30_q_le_10_matching_count"], 1)
        self.assertTrue(summary["candidate_grid_unique_match"])
        self.assertEqual(payload["matching_rational_candidates"], ["9/5"])
        self.assertTrue(summary["required_linear_gamma_strictly_increases"])
        self.assertTrue(summary["legacy_beta_tors_not_equal_any_required_gamma"])
        self.assertTrue(summary["damping_separation_inherited"])
        self.assertTrue(summary["exact_damping_calculus_inherited"])
        self.assertTrue(summary["finite_diagonal_completion_map_inherited"])
        self.assertTrue(summary["component_gap_damping_source_still_open"])
        self.assertFalse(summary["strict_beta_eta_source_exported"])
        self.assertFalse(summary["legacy_beta_tors_to_beta_eta_theorem_exported"])
        self.assertFalse(summary["full_bridge_theorem_exported"])

        self.assertTrue(all(payload["cross_checks"].values()))

    def test_proof_and_hard_limits(self):
        payload = self.payload
        proof = payload["proof_certificate"]
        self.assertIn("rg was used", proof["grep_step"])
        self.assertIn("S(1)-1=beta", proof["beta_step"])
        self.assertIn("log(d^(9/5))/log(d)=9/5", proof["eta_step"])
        self.assertIn("(S(d)-1)^5=d^9", proof["cover_step"])
        self.assertIn("exactly one matching candidate", proof["candidate_step"])
        self.assertIn("d^(4/5)", proof["legacy_linear_step"])
        self.assertIn("not a strict dynamical source", proof["theoretical_limit"])

        hard_limits = "\n".join(payload["hard_limits"])
        self.assertIn("No strict dynamical derivation of beta=1 or eta=9/5", hard_limits)
        self.assertIn("No beta_tors -> beta/eta theorem", hard_limits)
        self.assertIn("No beta_tors -> chi_11 theorem", hard_limits)
        self.assertIn("No legacy physical-role transfer", hard_limits)
        self.assertIn("No ToE closure", hard_limits)
        self.assertTrue(REPORT_MD.exists())


if __name__ == "__main__":
    unittest.main()
