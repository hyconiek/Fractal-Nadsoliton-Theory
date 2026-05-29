import json
import subprocess
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "fundamental_action_reconstruction" / "scratch" / "bridge_strict_alpha_hebbian_unit_orientation_bit_requirement_probe.py"
REPORT_JSON = ROOT / "fundamental_action_reconstruction" / "scratch" / "bridge_strict_alpha_hebbian_unit_orientation_bit_requirement_report.json"
REPORT_MD = ROOT / "fundamental_action_reconstruction" / "scratch" / "bridge_strict_alpha_hebbian_unit_orientation_bit_requirement_report.md"


class StrictAlphaHebbianUnitOrientationBitRequirementProbeTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        subprocess.run([sys.executable, str(SCRIPT)], check=True, cwd=ROOT)
        cls.payload = json.loads(REPORT_JSON.read_text(encoding="utf-8"))

    def test_result_kind_and_model(self):
        payload = self.payload
        self.assertEqual(
            payload["result_kind"],
            "SCRATCH_STRICT_ALPHA_HEBBIAN_UNIT_ORIENTATION_BIT_REQUIREMENT_PROBE__NOT_A_THEOREM",
        )
        self.assertEqual(payload["status"], "one-bit-unit-axis-record-required-for-singleton-d5-after-anti-nyquist-filtering")
        model = payload["finite_model"]
        self.assertEqual(model["ring"], "Z_12")
        self.assertEqual(model["support_count"], 792)
        self.assertEqual(model["histogram_class_count"], 35)
        self.assertEqual(model["automorphism_units"], [1, 5, 7, 11])
        self.assertEqual(model["target_q_power"], "256/243")
        self.assertEqual(model["target_eta"], "9/5")

    def test_character_table(self):
        rows = self.payload["axis_orbit"]["unit_action_character_table"]
        by_unit = {row["unit"]: row for row in rows}
        self.assertEqual(by_unit[1]["axis_swap_action"], "preserve")
        self.assertEqual(by_unit[11]["axis_swap_action"], "preserve")
        self.assertEqual(by_unit[5]["axis_swap_action"], "swap")
        self.assertEqual(by_unit[7]["axis_swap_action"], "swap")
        self.assertEqual([unit for unit, row in by_unit.items() if row["unit_orientation_character_chi_11_kernel"] == 1], [1, 11])
        self.assertEqual([unit for unit, row in by_unit.items() if row["unit_orientation_character_chi_11_kernel"] == -1], [5, 7])
        self.assertEqual(by_unit[5]["maps_k1_to"], "A5_k5_d5")
        self.assertEqual(by_unit[5]["maps_k5_to"], "A1_k1_contiguous")

    def test_boolean_selector_enumeration(self):
        enumeration = self.payload["boolean_selector_enumeration"]
        self.assertEqual(enumeration["selector_count"], 4)
        self.assertEqual(enumeration["full_aut_invariant_selector_names"], ["none", "both"])
        self.assertEqual(enumeration["stabilizer_1_11_invariant_selector_names"], ["none", "k1_only", "k5_d5_only", "both"])
        self.assertFalse(enumeration["full_aut_singleton_d5_exists"])
        self.assertTrue(enumeration["stabilizer_1_11_singleton_d5_exists"])

    def test_subgroup_requirement_certificate(self):
        cert = self.payload["subgroup_requirement_certificate"]
        self.assertEqual(cert["minimal_nontrivial_record"], "one binary unit-axis bit separating {1,11} from {5,7}")
        self.assertEqual(cert["stabilizer_kernel_for_d5_axis"], [1, 11])
        self.assertEqual(cert["swapping_coset"], [5, 7])
        self.assertEqual(cert["bits_required_on_two_axis_orbit"], 1)
        catalog = {tuple(row["subgroup"]): row for row in cert["subgroup_catalog"]}
        self.assertTrue(catalog[(1, 11)]["singleton_d5_invariant"])
        self.assertFalse(catalog[(1, 5, 7, 11)]["singleton_d5_invariant"])

    def test_interpretation_and_guardrails(self):
        interpretation = self.payload["selector_interpretation"]
        self.assertIn("one exact binary", interpretation["finite_gain"])
        self.assertIn("Full Aut(Z_12)", interpretation["negative_result"])
        self.assertIn("not derive", interpretation["honest_limit"])

        ontology = self.payload["ontology_guardrail"]
        self.assertIn("nadsoliton itself", ontology["allowed_reading"])
        self.assertIn("No separate informational layer", ontology["forbidden_reading"])

        hard_limits = "\n".join(self.payload["hard_limits"])
        self.assertIn("No identity K_legacy_ont == K_strict_gate", hard_limits)
        self.assertIn("No theorem derives the required one-bit", hard_limits)
        self.assertIn("No QW-2191 discharge", hard_limits)
        self.assertIn("No ToE closure", hard_limits)
        self.assertTrue(REPORT_MD.exists())


if __name__ == "__main__":
    unittest.main()
