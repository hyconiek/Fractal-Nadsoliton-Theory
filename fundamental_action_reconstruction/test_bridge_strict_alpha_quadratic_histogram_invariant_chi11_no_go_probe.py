import json
import subprocess
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "fundamental_action_reconstruction" / "scratch" / "bridge_strict_alpha_quadratic_histogram_invariant_chi11_no_go_probe.py"
REPORT_JSON = ROOT / "fundamental_action_reconstruction" / "scratch" / "bridge_strict_alpha_quadratic_histogram_invariant_chi11_no_go_report.json"
REPORT_MD = ROOT / "fundamental_action_reconstruction" / "scratch" / "bridge_strict_alpha_quadratic_histogram_invariant_chi11_no_go_report.md"


class StrictAlphaQuadraticHistogramInvariantChi11NoGoProbeTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        subprocess.run([sys.executable, str(SCRIPT)], check=True, cwd=ROOT)
        cls.payload = json.loads(REPORT_JSON.read_text(encoding="utf-8"))

    def test_result_kind_and_model(self):
        payload = self.payload
        self.assertEqual(
            payload["result_kind"],
            "SCRATCH_STRICT_ALPHA_QUADRATIC_HISTOGRAM_INVARIANT_CHI11_NO_GO_PROBE__NO_GO_NOT_A_THEOREM",
        )
        self.assertEqual(
            payload["status"],
            "degree2-full-aut-histogram-invariants-cannot-distinguish-a1-from-a5-or-export-chi11",
        )
        model = payload["finite_model"]
        self.assertEqual(model["ring"], "Z_12")
        self.assertEqual(model["active_count"], 5)
        self.assertEqual(model["automorphism_units"], [1, 5, 7, 11])
        self.assertEqual(model["branch_modes"], [1, 5, 7, 11])
        self.assertEqual(model["shell_labels"], ["d1", "d2", "d3", "d4", "d5", "d6"])
        self.assertEqual(model["monomial_count_degree_le_2"], 27)
        self.assertEqual(model["invariant_orbit_sum_count_degree_le_2"], 21)
        self.assertEqual(model["target_q_power"], "256/243")
        self.assertEqual(model["target_eta"], "9/5")

    def test_branch_rows_and_symbolic_certificate(self):
        branches = {row["name"]: row for row in self.payload["branch_rows"]}
        self.assertEqual(branches["A1_k1"]["distance_histogram_d1_to_d6"], [4, 3, 2, 1, 0, 0])
        self.assertEqual(branches["A5_k5"]["distance_histogram_d1_to_d6"], [0, 3, 2, 1, 4, 0])
        self.assertEqual(branches["A1_k1"]["required_chi11_value"], -1)
        self.assertEqual(branches["A5_k5"]["required_chi11_value"], 1)

        symbolic = self.payload["symbolic_certificate"]
        self.assertIn("d1<->d5 swap", symbolic["histogram_relation"])
        self.assertIn("P(A1)=P(A5)", symbolic["invariant_principle"])
        self.assertIn("6 linear + 21 quadratic", symbolic["degree2_enumeration"])
        self.assertIn("imports the missing shell-label bit", symbolic["anti_invariant_boundary"])

    def test_invariant_orbit_rows(self):
        rows = self.payload["invariant_orbit_rows"]
        self.assertEqual(len(rows), 21)
        self.assertTrue(all(not row["distinguishes_branch_pairs"] for row in rows))
        by_orbit = {tuple(row["orbit_monomials"]): row for row in rows}
        self.assertEqual(by_orbit[("d1", "d5")]["invariant_orbit_sum_value_A1_A11"], 4)
        self.assertEqual(by_orbit[("d1", "d5")]["invariant_orbit_sum_value_A5_A7"], 4)
        self.assertEqual(by_orbit[("d1*d1", "d5*d5")]["invariant_orbit_sum_value_A1_A11"], 16)
        self.assertEqual(by_orbit[("d1*d1", "d5*d5")]["invariant_orbit_sum_value_A5_A7"], 16)
        self.assertIn(("d1*d5",), by_orbit)

    def test_anti_invariant_rows_and_summary(self):
        summary = self.payload["quadratic_invariant_summary"]
        self.assertEqual(summary["invariant_orbit_sum_count"], 21)
        self.assertEqual(summary["invariant_orbit_sum_pair_distinguishing_count"], 0)
        self.assertEqual(summary["anti_invariant_orbit_difference_count"], 6)
        self.assertEqual(summary["anti_invariant_chi11_covariant_count"], 5)
        self.assertEqual(summary["allowed_strict_chi11_source_count"], 0)

        rows = {row["anti_orbit_difference"]: row for row in self.payload["anti_invariant_rows"]}
        self.assertEqual(rows["d5 - d1"]["value_A1_A11"], -4)
        self.assertEqual(rows["d5 - d1"]["value_A5_A7"], 4)
        self.assertTrue(rows["d5 - d1"]["numeric_chi11_covariant"])
        self.assertTrue(rows["d5*d5 - d1*d1"]["numeric_chi11_covariant"])
        self.assertFalse(rows["d5*d6 - d1*d6"]["numeric_chi11_covariant"])
        self.assertTrue(all(row["imports_d1_vs_d5_shell_orientation"] for row in rows.values()))

    def test_proof_interpretation_and_guardrails(self):
        proof = self.payload["exact_proof_certificate"]
        self.assertIn("degree<=2 histogram monomials", proof["finite_domain"])
        self.assertIn("21 invariant orbit sums", proof["invariant_no_go"])
        self.assertIn("P(h)=P(swap_d1_d5(h))", proof["symbolic_reason"])
        self.assertIn("Five anti-invariant", proof["anti_invariant_boundary"])
        self.assertIn("degree<=2 histogram invariants", proof["scope_limit"])

        interpretation = self.payload["interpretation"]
        self.assertIn("anti-invariant expressions", interpretation["honest_positive"])
        self.assertIn("No full-Aut invariant", interpretation["honest_negative"])
        self.assertIn("quadratic invariant orbit-sum basis", interpretation["relation_to_previous_probe"])

        ontology = self.payload["ontology_guardrail"]
        self.assertIn("nadsoliton itself", ontology["allowed_reading"])
        self.assertIn("No separate informational layer", ontology["forbidden_reading"])

        hard_limits = "\n".join(self.payload["hard_limits"])
        self.assertIn("No identity K_legacy_ont == K_strict_gate", hard_limits)
        self.assertIn("No theorem derives the chi_11-kernel", hard_limits)
        self.assertIn("degree<=2 histogram invariant no-go", hard_limits)
        self.assertIn("No QW-2191 discharge", hard_limits)
        self.assertIn("No ToE closure", hard_limits)
        self.assertTrue(REPORT_MD.exists())


if __name__ == "__main__":
    unittest.main()
