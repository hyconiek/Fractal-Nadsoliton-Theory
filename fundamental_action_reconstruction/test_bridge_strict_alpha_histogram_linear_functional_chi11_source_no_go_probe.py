import json
import subprocess
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "fundamental_action_reconstruction" / "scratch" / "bridge_strict_alpha_histogram_linear_functional_chi11_source_no_go_probe.py"
REPORT_JSON = ROOT / "fundamental_action_reconstruction" / "scratch" / "bridge_strict_alpha_histogram_linear_functional_chi11_source_no_go_report.json"
REPORT_MD = ROOT / "fundamental_action_reconstruction" / "scratch" / "bridge_strict_alpha_histogram_linear_functional_chi11_source_no_go_report.md"


class StrictAlphaHistogramLinearFunctionalChi11SourceNoGoProbeTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        subprocess.run([sys.executable, str(SCRIPT)], check=True, cwd=ROOT)
        cls.payload = json.loads(REPORT_JSON.read_text(encoding="utf-8"))

    def test_result_kind_and_model(self):
        payload = self.payload
        self.assertEqual(
            payload["result_kind"],
            "SCRATCH_STRICT_ALPHA_HISTOGRAM_LINEAR_FUNCTIONAL_CHI11_SOURCE_NO_GO_PROBE__NO_GO_NOT_A_THEOREM",
        )
        self.assertEqual(
            payload["status"],
            "bounded-shell-linear-histogram-functionals-export-chi11-only-by-importing-d1-vs-d5-label",
        )
        model = payload["finite_model"]
        self.assertEqual(model["ring"], "Z_12")
        self.assertEqual(model["active_count"], 5)
        self.assertEqual(model["automorphism_units"], [1, 5, 7, 11])
        self.assertEqual(model["branch_modes"], [1, 5, 7, 11])
        self.assertEqual(model["coefficient_domain"], [-1, 0, 1])
        self.assertEqual(model["weight_count"], 729)
        self.assertEqual(model["target_q_power"], "256/243")
        self.assertEqual(model["target_eta"], "9/5")

    def test_branch_rows(self):
        rows = {row["name"]: row for row in self.payload["branch_rows"]}
        self.assertEqual(rows["A1_k1"]["distance_histogram_d1_to_d6"], [4, 3, 2, 1, 0, 0])
        self.assertEqual(rows["A11_k11"]["distance_histogram_d1_to_d6"], [4, 3, 2, 1, 0, 0])
        self.assertEqual(rows["A5_k5"]["distance_histogram_d1_to_d6"], [0, 3, 2, 1, 4, 0])
        self.assertEqual(rows["A7_k7"]["distance_histogram_d1_to_d6"], [0, 3, 2, 1, 4, 0])
        self.assertEqual(rows["A1_k1"]["required_chi11_value"], -1)
        self.assertEqual(rows["A5_k5"]["required_chi11_value"], 1)

    def test_symbolic_certificate(self):
        symbolic = self.payload["symbolic_certificate"]
        self.assertIn("L(A1)=4*w1+3*w2+2*w3+w4", symbolic["score_A1"])
        self.assertIn("L(A5)=3*w2+2*w3+w4+4*w5", symbolic["score_A5"])
        self.assertIn("L(A5)-L(A1)=4*(w5-w1)", symbolic["pair_difference"])
        self.assertIn("L(A5)=-L(A1)", symbolic["chi11_condition"])
        self.assertIn("w1=w5", symbolic["full_aut_no_go"])

    def test_linear_functional_summary(self):
        summary = self.payload["linear_functional_summary"]
        self.assertEqual(summary["total_weight_count"], 729)
        self.assertEqual(summary["full_Aut_invariant_weight_count"], 243)
        self.assertEqual(summary["pair_distinguishing_weight_count"], 486)
        self.assertEqual(summary["full_Aut_pair_distinguishing_weight_count"], 0)
        self.assertEqual(summary["numeric_chi11_covariant_weight_count"], 54)
        self.assertEqual(summary["numeric_chi11_covariant_full_Aut_invariant_weight_count"], 0)
        self.assertEqual(summary["allowed_strict_chi11_source_weight_count"], 0)
        self.assertEqual(summary["chi11_covariant_imports_shell_label_count"], 54)
        self.assertEqual(summary["chi11_covariant_support_size_histogram"], {"2": 6, "3": 16, "4": 12, "5": 12, "6": 8})

    def test_minimal_chi11_examples_and_no_allowed_rows(self):
        examples = self.payload["minimal_chi11_covariant_examples"]
        self.assertEqual(len(examples), 12)
        first = examples[0]
        self.assertEqual(first["weights_d1_to_d6"], [-1, 0, 0, 0, 1, 0])
        self.assertEqual(first["score_A1_A11"], -4)
        self.assertEqual(first["score_A5_A7"], 4)
        self.assertEqual(first["score_sum_A1_plus_A5"], 0)
        self.assertTrue(first["numeric_chi11_covariant"])
        self.assertTrue(first["imports_d1_vs_d5_shell_label"])
        self.assertFalse(first["full_Aut_invariant_weight"])
        self.assertEqual(self.payload["allowed_strict_chi11_source_rows"], [])

    def test_proof_interpretation_and_guardrails(self):
        proof = self.payload["exact_proof_certificate"]
        self.assertIn("3^6=729", proof["finite_domain"])
        self.assertIn("L(A5)-L(A1)=4*(w5-w1)", proof["symbolic_no_go"])
        self.assertIn("allowed strict chi_11 source weights: 0", proof["enumerated_no_go"])
        self.assertIn("54 chi_11-covariant weights", proof["import_boundary"])
        self.assertIn("bounded shell-linear histogram", proof["scope_limit"])

        interpretation = self.payload["interpretation"]
        self.assertIn("many chi_11-covariant witnesses", interpretation["honest_positive"])
        self.assertIn("w1!=w5", interpretation["honest_negative"])
        self.assertIn("exhaustive bounded linear histogram audit", interpretation["relation_to_previous_probe"])

        ontology = self.payload["ontology_guardrail"]
        self.assertIn("nadsoliton itself", ontology["allowed_reading"])
        self.assertIn("No separate informational layer", ontology["forbidden_reading"])

        hard_limits = "\n".join(self.payload["hard_limits"])
        self.assertIn("No identity K_legacy_ont == K_strict_gate", hard_limits)
        self.assertIn("No theorem derives the chi_11-kernel", hard_limits)
        self.assertIn("bounded shell-linear histogram no-go", hard_limits)
        self.assertIn("No QW-2191 discharge", hard_limits)
        self.assertIn("No ToE closure", hard_limits)
        self.assertTrue(REPORT_MD.exists())


if __name__ == "__main__":
    unittest.main()
