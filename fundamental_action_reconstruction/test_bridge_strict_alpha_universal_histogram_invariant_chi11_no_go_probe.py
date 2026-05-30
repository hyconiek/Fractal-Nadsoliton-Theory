import json
import subprocess
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "fundamental_action_reconstruction" / "scratch" / "bridge_strict_alpha_universal_histogram_invariant_chi11_no_go_probe.py"
REPORT_JSON = ROOT / "fundamental_action_reconstruction" / "scratch" / "bridge_strict_alpha_universal_histogram_invariant_chi11_no_go_report.json"
REPORT_MD = ROOT / "fundamental_action_reconstruction" / "scratch" / "bridge_strict_alpha_universal_histogram_invariant_chi11_no_go_report.md"


class StrictAlphaUniversalHistogramInvariantChi11NoGoProbeTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        subprocess.run([sys.executable, str(SCRIPT)], check=True, cwd=ROOT)
        cls.payload = json.loads(REPORT_JSON.read_text(encoding="utf-8"))

    def test_result_kind_and_model(self):
        payload = self.payload
        self.assertEqual(
            payload["result_kind"],
            "SCRATCH_STRICT_ALPHA_UNIVERSAL_HISTOGRAM_INVARIANT_CHI11_NO_GO_PROBE__NO_GO_NOT_A_THEOREM",
        )
        self.assertEqual(
            payload["status"],
            "full-aut-invariant-histogram-only-sources-cannot-distinguish-a1-from-a5-or-export-chi11",
        )
        model = payload["finite_model"]
        self.assertEqual(model["ring"], "Z_12")
        self.assertEqual(model["active_count"], 5)
        self.assertEqual(model["automorphism_units"], [1, 5, 7, 11])
        self.assertEqual(model["support_count"], 792)
        self.assertEqual(model["histogram_class_count"], 35)
        self.assertEqual(model["full_Aut_histogram_orbit_count"], 24)
        self.assertEqual(model["target_q_power"], "256/243")
        self.assertEqual(model["target_eta"], "9/5")

    def test_a1_a5_orbit_certificate(self):
        cert = self.payload["a1_a5_histogram_orbit_certificate"]
        self.assertEqual(cert["A1_contiguous_histogram"], [4, 3, 2, 1, 0, 0])
        self.assertEqual(cert["A5_d5_histogram"], [0, 3, 2, 1, 4, 0])
        self.assertEqual(cert["swap_d1_d5_A1_histogram"], [0, 3, 2, 1, 4, 0])
        self.assertTrue(cert["same_full_Aut_histogram_orbit"])
        self.assertEqual(cert["A1_support_count"], 12)
        self.assertEqual(cert["A5_support_count"], 12)
        self.assertEqual(cert["orbit_support_count"], 24)

    def test_histogram_orbit_summary_and_rows(self):
        summary = self.payload["histogram_orbit_summary"]
        self.assertEqual(summary["histogram_class_count"], 35)
        self.assertEqual(summary["full_Aut_histogram_orbit_count"], 24)
        self.assertEqual(summary["fixed_histogram_orbit_count"], 13)
        self.assertEqual(summary["two_histogram_orbit_count"], 11)
        self.assertEqual(summary["invariant_boolean_classifier_count_power"], 24)
        self.assertEqual(summary["invariant_boolean_classifier_total"], 16777216)
        self.assertEqual(summary["singleton_A5_not_A1_invariant_classifier_count"], 0)
        self.assertEqual(summary["both_or_neither_A1_A5_invariant_classifier_count"], 16777216)

        rows = self.payload["histogram_orbit_rows"]
        self.assertEqual(len(rows), 24)
        a1_a5_rows = [row for row in rows if row["is_A1_A5_pair_orbit"]]
        self.assertEqual(len(a1_a5_rows), 1)
        row = a1_a5_rows[0]
        self.assertEqual(row["orbit_size"], 2)
        self.assertEqual(row["support_count_total"], 24)
        self.assertTrue(row["contains_A1_contiguous_histogram"])
        self.assertTrue(row["contains_A5_d5_histogram"])

    def test_proof_interpretation_and_guardrails(self):
        proof = self.payload["exact_proof_certificate"]
        self.assertIn("C(12,5)=792", proof["finite_domain"])
        self.assertIn("d1<->d5 swap", proof["histogram_swap_fact"])
        self.assertIn("F(A1)=F(A5)", proof["universal_invariant_no_go"])
        self.assertIn("35 histogram classes", proof["classifier_count"])
        self.assertIn("singleton A5-vs-A1 classification count is 0", proof["classifier_count"])
        self.assertIn("imports the missing shell-label orientation", proof["anti_invariant_boundary"])
        self.assertIn("histogram-only full-Aut invariant sources", proof["scope_limit"])

        interpretation = self.payload["interpretation"]
        self.assertIn("every possible full-Aut invariant function", interpretation["honest_positive"])
        self.assertIn("non-histogram strict sources open", interpretation["honest_negative"])
        self.assertIn("subsumes the linear and quadratic", interpretation["relation_to_previous_probe"])

        ontology = self.payload["ontology_guardrail"]
        self.assertIn("nadsoliton itself", ontology["allowed_reading"])
        self.assertIn("No separate informational layer", ontology["forbidden_reading"])

        hard_limits = "\n".join(self.payload["hard_limits"])
        self.assertIn("No identity K_legacy_ont == K_strict_gate", hard_limits)
        self.assertIn("No theorem derives the chi_11-kernel", hard_limits)
        self.assertIn("universal histogram-only full-Aut invariant no-go", hard_limits)
        self.assertIn("No QW-2191 discharge", hard_limits)
        self.assertIn("No ToE closure", hard_limits)
        self.assertTrue(REPORT_MD.exists())


if __name__ == "__main__":
    unittest.main()
