import json
import subprocess
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "fundamental_action_reconstruction" / "scratch" / "bridge_strict_alpha_d12_chi11_max_shell_imbalance_selector_certificate_probe.py"
REPORT_JSON = ROOT / "fundamental_action_reconstruction" / "scratch" / "bridge_strict_alpha_d12_chi11_max_shell_imbalance_selector_certificate_report.json"
REPORT_MD = ROOT / "fundamental_action_reconstruction" / "scratch" / "bridge_strict_alpha_d12_chi11_max_shell_imbalance_selector_certificate_report.md"


class StrictAlphaD12Chi11MaxShellImbalanceSelectorCertificateProbeTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        subprocess.run([sys.executable, str(SCRIPT)], check=True, cwd=ROOT)
        cls.payload = json.loads(REPORT_JSON.read_text(encoding="utf-8"))

    def test_result_kind_and_finite_model(self):
        payload = self.payload
        self.assertEqual(
            payload["result_kind"],
            "SCRATCH_STRICT_ALPHA_D12_CHI11_MAX_SHELL_IMBALANCE_SELECTOR_CERTIFICATE__CONDITIONAL_NOT_A_THEOREM",
        )
        self.assertEqual(
            payload["status"],
            "branch-generator-uniquely-maximizes-shell-labelled-d1-d5-imbalance-but-imports-axis",
        )
        model = payload["finite_model"]
        self.assertEqual(model["ring"], "Z_12")
        self.assertEqual(model["active_count"], 5)
        self.assertEqual(model["support_count"], 792)
        self.assertEqual(model["dihedral_units"], [1, 11])
        self.assertEqual(model["automorphism_units"], [1, 5, 7, 11])
        self.assertEqual(model["d12_orbit_count"], 38)
        self.assertEqual(model["unit5_two_cycle_count"], 13)
        self.assertEqual(model["target_q_power"], "256/243")
        self.assertEqual(model["target_eta"], "9/5")

    def test_selector_summary_and_cycle_scores(self):
        summary = self.payload["selector_summary"]
        self.assertEqual(summary["candidate_score"], "abs(h_d5-h_d1) on each D12/unit5 two-cycle")
        self.assertTrue(summary["score_requires_shell_label"])
        self.assertTrue(summary["score_requires_reduced_d12_quotient"])
        self.assertEqual(summary["max_absolute_shell_imbalance"], 4)
        self.assertEqual(summary["maximizer_count"], 1)
        self.assertTrue(summary["unique_maximizer_is_branch_cycle"])
        self.assertEqual(summary["absolute_shell_imbalance_distribution"], {"0": 2, "1": 6, "2": 4, "4": 1})
        self.assertFalse(summary["full_aut_allowed_strict_source"])

        rows = self.payload["cycle_score_rows"]
        self.assertEqual(len(rows), 13)
        self.assertTrue(all(row["imbalance_flips_under_unit5"] for row in rows))
        self.assertEqual(sum(1 for row in rows if row["absolute_shell_imbalance"] == 4), 1)
        self.assertEqual(sum(1 for row in rows if row["is_branch_cycle"]), 1)
        self.assertEqual(
            sorted(row["absolute_shell_imbalance"] for row in rows),
            [0, 0, 1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 4],
        )

    def test_branch_maximizer_and_branch_rows(self):
        cert = self.payload["branch_maximizer_certificate"]
        self.assertEqual(cert["cycle_index"], 0)
        self.assertEqual(cert["orbit_pair"], [0, 37])
        self.assertEqual(cert["low_histogram_d1_to_d6"], [4, 3, 2, 1, 0, 0])
        self.assertEqual(cert["high_histogram_d1_to_d6"], [0, 3, 2, 1, 4, 0])
        self.assertEqual(cert["low_d5_minus_d1"], -4)
        self.assertEqual(cert["high_d5_minus_d1"], 4)
        self.assertEqual(cert["absolute_shell_imbalance"], 4)
        self.assertEqual(cert["low_gap_necklace"], [1, 1, 1, 1, 8])
        self.assertEqual(cert["high_gap_necklace"], [2, 2, 3, 2, 3])
        self.assertTrue(cert["is_branch_cycle"])

        branches = {row["name"]: row for row in self.payload["branch_rows"]}
        self.assertEqual(branches["A1_k1"]["d5_minus_d1"], -4)
        self.assertEqual(branches["A11_k11"]["d5_minus_d1"], -4)
        self.assertEqual(branches["A5_k5"]["d5_minus_d1"], 4)
        self.assertEqual(branches["A7_k7"]["d5_minus_d1"], 4)
        self.assertEqual(branches["A1_k1"]["required_chi11_value"], -1)
        self.assertEqual(branches["A5_k5"]["required_chi11_value"], 1)

    def test_proof_interpretation_and_guardrails(self):
        proof = self.payload["exact_proof_certificate"]
        self.assertIn("C(12,5)=792", proof["finite_domain"])
        self.assertIn("I(unit5*O)=-I(O)", proof["covariance_fact"])
        self.assertIn("max |I| is 4", proof["maximality_certificate"])
        self.assertIn("attained by 1 cycle", proof["maximality_certificate"])
        self.assertIn("unique maximizer", proof["branch_selection"])
        self.assertIn("labelled d1/d5 shell axis", proof["import_boundary"])
        self.assertIn("not an invariant source", proof["full_aut_intersection"])

        interpretation = self.payload["interpretation"]
        self.assertIn("unique maximum shell-imbalance", interpretation["honest_positive"])
        self.assertIn("imports the d1/d5 shell axis", interpretation["honest_negative"])
        self.assertIn("selecting one generator", interpretation["relation_to_previous_probe"])

        ontology = self.payload["ontology_guardrail"]
        self.assertIn("nadsoliton itself", ontology["allowed_reading"])
        self.assertIn("No separate informational layer", ontology["forbidden_reading"])

        hard_limits = "\n".join(self.payload["hard_limits"])
        self.assertIn("No identity K_legacy_ont == K_strict_gate", hard_limits)
        self.assertIn("No theorem derives chi_11", hard_limits)
        self.assertIn("No QW-2191 discharge", hard_limits)
        self.assertIn("No ToE closure", hard_limits)
        self.assertIn("not a full-Aut strict-source theorem", hard_limits)
        self.assertTrue(REPORT_MD.exists())


if __name__ == "__main__":
    unittest.main()
