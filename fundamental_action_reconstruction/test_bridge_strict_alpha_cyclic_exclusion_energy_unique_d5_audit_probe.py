import json
import subprocess
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "fundamental_action_reconstruction" / "scratch" / "bridge_strict_alpha_cyclic_exclusion_energy_unique_d5_audit_probe.py"
REPORT_JSON = ROOT / "fundamental_action_reconstruction" / "scratch" / "bridge_strict_alpha_cyclic_exclusion_energy_unique_d5_audit_report.json"
REPORT_MD = ROOT / "fundamental_action_reconstruction" / "scratch" / "bridge_strict_alpha_cyclic_exclusion_energy_unique_d5_audit_report.md"


class StrictAlphaCyclicExclusionEnergyUniqueD5AuditProbeTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        subprocess.run([sys.executable, str(SCRIPT)], check=True, cwd=ROOT)
        cls.payload = json.loads(REPORT_JSON.read_text(encoding="utf-8"))

    def test_result_kind_and_model(self):
        payload = self.payload
        self.assertEqual(
            payload["result_kind"],
            "SCRATCH_STRICT_ALPHA_CYCLIC_EXCLUSION_ENERGY_UNIQUE_D5_AUDIT_PROBE__CONDITIONAL_NOT_A_THEOREM",
        )
        self.assertEqual(
            payload["status"],
            "d1_plus_d6_exclusion_energy-uniquely-selects-A5-d5-orbit-conditionally-source-not-derived",
        )
        model = payload["finite_model"]
        self.assertEqual(model["ring"], "Z_12")
        self.assertEqual(model["support_count"], 792)
        self.assertEqual(model["histogram_class_count"], 35)
        self.assertEqual(model["target_q_power"], "256/243")
        self.assertEqual(model["target_eta"], "9/5")

    def test_pairwise_and_exclusion_energy_winners(self):
        pairwise = self.payload["pairwise_A1_A5_audit"]
        self.assertEqual(pairwise["A1_histogram_d1_to_d6"], [4, 3, 2, 1, 0, 0])
        self.assertEqual(pairwise["A1_exclusion_energy_d1_plus_d6"], 4)
        self.assertEqual(pairwise["A5_histogram_d1_to_d6"], [0, 3, 2, 1, 4, 0])
        self.assertEqual(pairwise["A5_exclusion_energy_d1_plus_d6"], 0)
        self.assertTrue(pairwise["exclusion_energy_prefers_A5_over_A1_pairwise"])

        energy = self.payload["cyclic_exclusion_energy_audit"]
        self.assertEqual(energy["energy_definition"], "E_excl = d1_nearest + d6_antipodal")
        self.assertEqual(energy["minimum_energy"], 0)
        self.assertEqual(energy["winner_support_count"], 12)
        self.assertEqual(energy["winner_histogram_class_count"], 1)
        self.assertEqual(energy["winner_dihedral_orbit_count"], 1)
        self.assertEqual(energy["winner_histograms"], [{"distance_histogram_d1_to_d6": [0, 3, 2, 1, 4, 0], "support_count": 12}])
        self.assertTrue(energy["selects_A5_d5_orbit"])
        self.assertFalse(energy["selects_A1_contiguous_orbit"])
        self.assertEqual(energy["winner_orbits"][0]["exclusion_energy_d1_plus_d6"], 0)
        self.assertTrue(energy["winner_orbits"][0]["is_A5_d5_orbit"])

    def test_energy_distribution_and_d1_refinement(self):
        energy_distribution = {
            row["energy"]: row["support_count"]
            for row in self.payload["cyclic_exclusion_energy_audit"]["energy_distribution"]
        }
        self.assertEqual(energy_distribution, {0: 12, 1: 60, 2: 228, 3: 348, 4: 120, 5: 24})

        d1_refinement = self.payload["d1_zero_refinement_audit"]
        self.assertEqual(d1_refinement["d1_zero_support_count"], 36)
        self.assertEqual(d1_refinement["d1_zero_histogram_class_count"], 3)
        self.assertEqual(d1_refinement["d1_zero_dihedral_orbit_count"], 3)
        rows = d1_refinement["d1_zero_histogram_rows"]
        self.assertEqual([row["d6_antipodal_count"] for row in rows], [0, 1, 2])
        self.assertEqual(sum(row["survives_d1_plus_d6_zero"] for row in rows), 1)
        self.assertEqual(sum(row["is_A5_d5_histogram"] for row in rows), 1)

    def test_binary_shell_mask_classification(self):
        mask_audit = self.payload["binary_shell_mask_audit"]
        self.assertEqual(mask_audit["binary_masks_tested"], 63)
        self.assertEqual(mask_audit["A5_selecting_binary_mask_count"], 7)
        self.assertTrue(mask_audit["contains_d1_plus_d6_mask"])
        minimal_masks = {
            tuple(row["penalized_shell_indices_1_based"])
            for row in mask_audit["inclusion_minimal_A5_selecting_masks"]
        }
        self.assertEqual(minimal_masks, {(1, 4), (1, 6)})
        d1_d6 = next(
            row for row in mask_audit["A5_selecting_binary_masks"]
            if row["penalized_shell_indices_1_based"] == [1, 6]
        )
        self.assertEqual(d1_d6["penalized_shell_names"], ["d1_nearest", "d6_antipodal"])
        self.assertEqual(d1_d6["minimum_energy"], 0)
        self.assertEqual(d1_d6["winner_dihedral_orbit_count"], 1)
        self.assertIn("not a derivation", mask_audit["honest_status"])

    def test_proof_interpretation_and_guardrails(self):
        proof = self.payload["exact_proof_certificate"]
        self.assertIn("C(12,5)=792", proof["enumeration_domain"])
        self.assertIn("d1(S)+d6(S)", proof["energy_definition"])
        self.assertIn("one dihedral orbit", proof["unique_orbit_fact"])
        self.assertIn("36 supports", proof["d1_refinement"])
        self.assertIn("63 nonempty binary shell penalties", proof["binary_mask_classification"])
        self.assertIn("not a derivation", proof["missing_source"])

        interpretation = self.payload["interpretation"]
        self.assertIn("uniquely selects", interpretation["direct_result"])
        self.assertIn("adding antipodal exclusion", interpretation["improvement_over_d1_only"])
        self.assertIn("cleaner shell-exclusion witness", interpretation["relation_to_previous_tie_break"])
        self.assertIn("conditional/non-strict", interpretation["honest_limit"])

        ontology = self.payload["ontology_guardrail"]
        self.assertIn("nadsoliton itself", ontology["allowed_reading"])
        self.assertIn("No separate informational layer", ontology["forbidden_reading"])

        hard_limits = "\n".join(self.payload["hard_limits"])
        self.assertIn("No identity K_legacy_ont == K_strict_gate", hard_limits)
        self.assertIn("No theorem derives cyclic exclusion energy", hard_limits)
        self.assertIn("conditional finite witness", hard_limits)
        self.assertIn("No QW-2191 discharge", hard_limits)
        self.assertIn("No ToE closure", hard_limits)
        self.assertTrue(REPORT_MD.exists())


if __name__ == "__main__":
    unittest.main()
