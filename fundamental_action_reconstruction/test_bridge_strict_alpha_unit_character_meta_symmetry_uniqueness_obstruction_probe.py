import json
import subprocess
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "fundamental_action_reconstruction" / "scratch" / "bridge_strict_alpha_unit_character_meta_symmetry_uniqueness_obstruction_probe.py"
REPORT_JSON = ROOT / "fundamental_action_reconstruction" / "scratch" / "bridge_strict_alpha_unit_character_meta_symmetry_uniqueness_obstruction_report.json"
REPORT_MD = ROOT / "fundamental_action_reconstruction" / "scratch" / "bridge_strict_alpha_unit_character_meta_symmetry_uniqueness_obstruction_report.md"


class StrictAlphaUnitCharacterMetaSymmetryUniquenessObstructionProbeTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        subprocess.run([sys.executable, str(SCRIPT)], check=True, cwd=ROOT)
        cls.payload = json.loads(REPORT_JSON.read_text(encoding="utf-8"))

    def test_result_kind_and_model(self):
        payload = self.payload
        self.assertEqual(
            payload["result_kind"],
            "SCRATCH_STRICT_ALPHA_UNIT_CHARACTER_META_SYMMETRY_UNIQUENESS_OBSTRUCTION_PROBE__NOT_A_THEOREM",
        )
        self.assertEqual(payload["status"], "abstract-nontrivial-character-class-does-not-uniquely-select-chi_11")
        model = payload["finite_model"]
        self.assertEqual(model["ring"], "Z_12")
        self.assertEqual(model["support_count"], 792)
        self.assertEqual(model["histogram_class_count"], 35)
        self.assertEqual(model["automorphism_units"], [1, 5, 7, 11])
        self.assertEqual(model["target_q_power"], "256/243")
        self.assertEqual(model["target_eta"], "9/5")

    def test_meta_automorphism_and_character_orbit(self):
        meta = self.payload["meta_automorphism_audit"]
        self.assertEqual(meta["meta_automorphism_count"], 6)
        self.assertTrue(all(row["unit_mapping"]["1"] == 1 for row in meta["rows"]))

        orbit_audit = self.payload["character_meta_orbit_audit"]
        self.assertEqual(orbit_audit["character_count"], 4)
        required = orbit_audit["required_chi_11_meta_orbit"]
        self.assertEqual(
            required["meta_orbit"],
            ["chi_11_required_d5_unit_axis", "chi_5_kernel", "chi_7_kernel"],
        )
        self.assertEqual(required["meta_orbit_size"], 3)
        self.assertEqual(orbit_audit["required_chi_11_stabilizer_size_in_meta_group"], 2)
        self.assertTrue(orbit_audit["nontrivial_characters_form_single_meta_orbit"])
        image_names = {row["required_chi_11_maps_to"] for row in orbit_audit["meta_action_on_required_chi_11"]}
        self.assertEqual(image_names, {"chi_11_required_d5_unit_axis", "chi_5_kernel", "chi_7_kernel"})

    def test_selection_rule_audit(self):
        by_rule = {row["rule"]: row for row in self.payload["selection_rule_audit"]}
        any_nontrivial = by_rule["choose_any_nontrivial_character"]
        self.assertTrue(any_nontrivial["meta_invariant"])
        self.assertFalse(any_nontrivial["selects_unique_required_chi_11"])
        self.assertEqual(len(any_nontrivial["selected_set"]), 3)

        kernel_size = by_rule["choose_character_with_kernel_size_2"]
        self.assertTrue(kernel_size["meta_invariant"])
        self.assertFalse(kernel_size["selects_unique_required_chi_11"])

        chi_11 = by_rule["choose_chi_11_by_kernel_{1,11}"]
        self.assertFalse(chi_11["meta_invariant"])
        self.assertTrue(chi_11["selects_unique_required_chi_11"])
        self.assertIn("unit-label/geometry", chi_11["honest_status"])

    def test_proof_interpretation_and_guardrails(self):
        proof = self.payload["exact_proof_certificate"]
        self.assertIn("6 elements", proof["meta_group"])
        self.assertIn("permutes the three nontrivial", proof["dual_action"])
        self.assertIn("meta-orbit of size 3", proof["orbit_fact"])
        self.assertIn("singleton chi_11 is not meta-invariant", proof["uniqueness_obstruction"])
        self.assertIn("{1,11} d5 stabilizer", proof["required_extra_datum"])

        interpretation = self.payload["interpretation"]
        self.assertIn("three-character ambiguity", interpretation["direct_result"])
        self.assertIn("No meta-invariant rule", interpretation["negative_result"])
        self.assertIn("identifies 11", interpretation["conditional_positive_result"])
        self.assertIn("does not derive", interpretation["honest_limit"])

        ontology = self.payload["ontology_guardrail"]
        self.assertIn("nadsoliton itself", ontology["allowed_reading"])
        self.assertIn("No separate informational layer", ontology["forbidden_reading"])

        hard_limits = "\n".join(self.payload["hard_limits"])
        self.assertIn("No identity K_legacy_ont == K_strict_gate", hard_limits)
        self.assertIn("No theorem derives chi_11", hard_limits)
        self.assertIn("does not discharge QW-2191", hard_limits)
        self.assertIn("No QW-2191 discharge", hard_limits)
        self.assertIn("No ToE closure", hard_limits)
        self.assertTrue(REPORT_MD.exists())


if __name__ == "__main__":
    unittest.main()
