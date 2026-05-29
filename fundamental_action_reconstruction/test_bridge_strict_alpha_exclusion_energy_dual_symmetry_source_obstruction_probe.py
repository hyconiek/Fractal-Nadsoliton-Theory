import json
import subprocess
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "fundamental_action_reconstruction" / "scratch" / "bridge_strict_alpha_exclusion_energy_dual_symmetry_source_obstruction_probe.py"
REPORT_JSON = ROOT / "fundamental_action_reconstruction" / "scratch" / "bridge_strict_alpha_exclusion_energy_dual_symmetry_source_obstruction_report.json"
REPORT_MD = ROOT / "fundamental_action_reconstruction" / "scratch" / "bridge_strict_alpha_exclusion_energy_dual_symmetry_source_obstruction_report.md"


class StrictAlphaExclusionEnergyDualSymmetrySourceObstructionProbeTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        subprocess.run([sys.executable, str(SCRIPT)], check=True, cwd=ROOT)
        cls.payload = json.loads(REPORT_JSON.read_text(encoding="utf-8"))

    def test_result_kind_and_model(self):
        payload = self.payload
        self.assertEqual(
            payload["result_kind"],
            "SCRATCH_STRICT_ALPHA_EXCLUSION_ENERGY_DUAL_SYMMETRY_SOURCE_OBSTRUCTION_PROBE__NOT_A_THEOREM",
        )
        self.assertEqual(payload["status"], "d1_plus_d6-selector-success-imports-chi_11-shell-label-premise")
        model = payload["finite_model"]
        self.assertEqual(model["ring"], "Z_12")
        self.assertEqual(model["support_count"], 792)
        self.assertEqual(model["histogram_class_count"], 35)
        self.assertEqual(model["automorphism_units"], [1, 5, 7, 11])
        self.assertEqual(model["target_q_power"], "256/243")
        self.assertEqual(model["target_eta"], "9/5")
        self.assertEqual(model["A1_histogram_d1_to_d6"], [4, 3, 2, 1, 0, 0])
        self.assertEqual(model["A5_histogram_d1_to_d6"], [0, 3, 2, 1, 4, 0])

    def test_shell_action_and_stabilizer(self):
        action = self.payload["shell_action_audit"]
        self.assertEqual(action["E_A5_formula"], "d1_nearest + d6_antipodal")
        self.assertEqual(action["stabilizer_units"], [1, 11])
        self.assertTrue(action["stabilizer_is_chi_11_kernel"])
        self.assertEqual(action["Aut_orbit_size"], 2)
        orbit_names = {row["name"] for row in action["Aut_orbit_weights"]}
        self.assertEqual(orbit_names, {"E_A5_candidate_d1_plus_d6", "E_A1_conjugate_d5_plus_d6"})
        rows_by_unit = {row["unit"]: row for row in action["rows"]}
        self.assertEqual(rows_by_unit[1]["shell_permutation_source_to_target"], [1, 2, 3, 4, 5, 6])
        self.assertEqual(rows_by_unit[11]["shell_permutation_source_to_target"], [1, 2, 3, 4, 5, 6])
        self.assertEqual(rows_by_unit[5]["shell_permutation_source_to_target"], [5, 2, 3, 4, 1, 6])
        self.assertEqual(rows_by_unit[7]["shell_permutation_source_to_target"], [5, 2, 3, 4, 1, 6])
        self.assertTrue(rows_by_unit[1]["preserves_E_A5"])
        self.assertFalse(rows_by_unit[5]["preserves_E_A5"])

    def test_aut_conjugate_selector_outcomes(self):
        by_name = {row["name"]: row for row in self.payload["aut_conjugate_selector_audit"]}
        a5 = by_name["E_A5_candidate_d1_plus_d6"]
        self.assertEqual(a5["minimum_score"], 0)
        self.assertEqual(a5["winner_support_count"], 12)
        self.assertEqual(a5["winner_dihedral_orbit_count"], 1)
        self.assertTrue(a5["selects_A5_d5_orbit"])
        self.assertFalse(a5["selects_A1_contiguous_orbit"])
        self.assertEqual(a5["winner_histogram_rows"], [{"distance_histogram_d1_to_d6": [0, 3, 2, 1, 4, 0], "support_count": 12}])

        a1 = by_name["E_A1_conjugate_d5_plus_d6"]
        self.assertEqual(a1["minimum_score"], 0)
        self.assertEqual(a1["winner_support_count"], 12)
        self.assertEqual(a1["winner_dihedral_orbit_count"], 1)
        self.assertTrue(a1["selects_A1_contiguous_orbit"])
        self.assertFalse(a1["selects_A5_d5_orbit"])
        self.assertEqual(a1["winner_histogram_rows"], [{"distance_histogram_d1_to_d6": [4, 3, 2, 1, 0, 0], "support_count": 12}])

    def test_proof_interpretation_and_guardrails(self):
        proof = self.payload["exact_proof_certificate"]
        self.assertIn("d1 to d5", proof["shell_action"])
        self.assertIn("{1,11}", proof["stabilizer"])
        self.assertIn("two weights", proof["orbit_fact"])
        self.assertIn("d5+d6 selects A1", proof["selector_conjugates"])
        self.assertIn("unit-label/chi_11 premise", proof["source_obstruction"])

        interpretation = self.payload["interpretation"]
        self.assertIn("not full-Aut invariant", interpretation["direct_result"])
        self.assertIn("selects A1/contiguous", interpretation["conjugate_warning"])
        self.assertIn("conditional", interpretation["honest_limit"])
        self.assertIn("symmetry cost", interpretation["research_value"])

        ontology = self.payload["ontology_guardrail"]
        self.assertIn("nadsoliton itself", ontology["allowed_reading"])
        self.assertIn("No separate informational layer", ontology["forbidden_reading"])

        hard_limits = "\n".join(self.payload["hard_limits"])
        self.assertIn("No identity K_legacy_ont == K_strict_gate", hard_limits)
        self.assertIn("No theorem derives the d1-vs-d5 shell label", hard_limits)
        self.assertIn("chi_11-kernel-stabilized conditional data", hard_limits)
        self.assertIn("No QW-2191 discharge", hard_limits)
        self.assertIn("No ToE closure", hard_limits)
        self.assertTrue(REPORT_MD.exists())


if __name__ == "__main__":
    unittest.main()
