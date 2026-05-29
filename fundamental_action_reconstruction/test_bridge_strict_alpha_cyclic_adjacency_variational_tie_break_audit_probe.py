import json
import subprocess
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "fundamental_action_reconstruction" / "scratch" / "bridge_strict_alpha_cyclic_adjacency_variational_tie_break_audit_probe.py"
REPORT_JSON = ROOT / "fundamental_action_reconstruction" / "scratch" / "bridge_strict_alpha_cyclic_adjacency_variational_tie_break_audit_report.json"
REPORT_MD = ROOT / "fundamental_action_reconstruction" / "scratch" / "bridge_strict_alpha_cyclic_adjacency_variational_tie_break_audit_report.md"


class StrictAlphaCyclicAdjacencyVariationalTieBreakAuditProbeTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        subprocess.run([sys.executable, str(SCRIPT)], check=True, cwd=ROOT)
        cls.payload = json.loads(REPORT_JSON.read_text(encoding="utf-8"))

    def test_result_kind_and_model(self):
        payload = self.payload
        self.assertEqual(
            payload["result_kind"],
            "SCRATCH_STRICT_ALPHA_CYCLIC_ADJACENCY_VARIATIONAL_TIE_BREAK_AUDIT_PROBE__CONDITIONAL_NOT_A_THEOREM",
        )
        self.assertEqual(
            payload["status"],
            "nearest-neighbor-adjacency-alone-not-global-d5-selector-lexicographic-d1-d5-rule-conditionally-selects-A5",
        )
        model = payload["finite_model"]
        self.assertEqual(model["ring"], "Z_12")
        self.assertEqual(model["support_count"], 792)
        self.assertEqual(model["histogram_class_count"], 35)
        self.assertEqual(model["target_q_power"], "256/243")
        self.assertEqual(model["target_eta"], "9/5")

    def test_pairwise_and_d1_distribution(self):
        pairwise = self.payload["pairwise_A1_A5_audit"]
        self.assertEqual(pairwise["A1_histogram_d1_to_d6"], [4, 3, 2, 1, 0, 0])
        self.assertEqual(pairwise["A5_histogram_d1_to_d6"], [0, 3, 2, 1, 4, 0])
        self.assertTrue(pairwise["minimize_d1_prefers_A5_over_A1_pairwise"])
        self.assertIn("weaker than global", pairwise["global_warning"])

        distribution = {row["d1_nearest_neighbor_count"]: row for row in self.payload["d1_distribution"]}
        self.assertEqual(distribution[0]["support_count"], 36)
        self.assertEqual(distribution[0]["histogram_class_count"], 3)
        self.assertEqual(distribution[0]["dihedral_orbit_count"], 3)
        self.assertEqual(distribution[4]["support_count"], 12)
        self.assertEqual(distribution[4]["dihedral_orbit_count"], 1)

    def test_adjacency_minimizer_no_go(self):
        minimizer = self.payload["adjacency_minimizer_audit"]
        self.assertEqual(minimizer["minimum_d1"], 0)
        self.assertEqual(minimizer["minimum_d1_support_count"], 36)
        self.assertEqual(minimizer["minimum_d1_histogram_class_count"], 3)
        self.assertEqual(minimizer["minimum_d1_dihedral_orbit_count"], 3)
        self.assertFalse(minimizer["nearest_neighbor_minimization_uniquely_selects_A5"])
        histogram_rows = minimizer["minimum_d1_histogram_rows"]
        self.assertEqual(len(histogram_rows), 3)
        self.assertEqual(sum(1 for row in histogram_rows if row["is_A5_d5_histogram"]), 1)
        self.assertIn(
            [0, 3, 2, 1, 4, 0],
            [row["distance_histogram_d1_to_d6"] for row in histogram_rows],
        )

    def test_lexicographic_tie_break_conditional_positive(self):
        tie_break = self.payload["lexicographic_tie_break_audit"]
        self.assertEqual(tie_break["rule"], "lexicographic_minimize_d1_then_maximize_d5")
        self.assertEqual(tie_break["best_d1"], 0)
        self.assertEqual(tie_break["best_d5"], 4)
        self.assertEqual(tie_break["winner_support_count"], 12)
        self.assertEqual(tie_break["winner_dihedral_orbit_count"], 1)
        self.assertEqual(tie_break["winner_histograms"], [[0, 3, 2, 1, 4, 0]])
        self.assertTrue(tie_break["selects_A5_d5_orbit"])
        self.assertIn("conditional", tie_break["honest_status"])
        self.assertTrue(tie_break["winner_orbits"][0]["is_A5_d5_orbit"])
        self.assertFalse(tie_break["winner_orbits"][0]["is_A1_contiguous_orbit"])

    def test_proof_interpretation_and_guardrails(self):
        proof = self.payload["exact_proof_certificate"]
        self.assertIn("C(12,5)=792", proof["enumeration_domain"])
        self.assertIn("d1", proof["adjacency_energy"])
        self.assertIn("3 dihedral orbits", proof["no_go"])
        self.assertIn("A5/d5 orbit", proof["conditional_positive"])
        self.assertIn("extra variational structure", proof["missing_source"])

        interpretation = self.payload["interpretation"]
        self.assertIn("does not globally select", interpretation["direct_result"])
        self.assertIn("A5 has d1=0", interpretation["pairwise_result"])
        self.assertIn("three dihedral-orbit minimizers", interpretation["global_negative_result"])
        self.assertIn("maximize d5", interpretation["conditional_positive_result"])
        self.assertIn("not derived from strict core", interpretation["honest_limit"])

        ontology = self.payload["ontology_guardrail"]
        self.assertIn("nadsoliton itself", ontology["allowed_reading"])
        self.assertIn("No separate informational layer", ontology["forbidden_reading"])

        hard_limits = "\n".join(self.payload["hard_limits"])
        self.assertIn("No identity K_legacy_ont == K_strict_gate", hard_limits)
        self.assertIn("No theorem derives cyclic adjacency", hard_limits)
        self.assertIn("conditional variational tie-break", hard_limits)
        self.assertIn("No QW-2191 discharge", hard_limits)
        self.assertIn("No ToE closure", hard_limits)
        self.assertTrue(REPORT_MD.exists())


if __name__ == "__main__":
    unittest.main()
