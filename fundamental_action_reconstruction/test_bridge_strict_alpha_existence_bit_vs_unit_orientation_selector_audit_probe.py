import json
import subprocess
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "fundamental_action_reconstruction" / "scratch" / "bridge_strict_alpha_existence_bit_vs_unit_orientation_selector_audit_probe.py"
REPORT_JSON = ROOT / "fundamental_action_reconstruction" / "scratch" / "bridge_strict_alpha_existence_bit_vs_unit_orientation_selector_audit_report.json"
REPORT_MD = ROOT / "fundamental_action_reconstruction" / "scratch" / "bridge_strict_alpha_existence_bit_vs_unit_orientation_selector_audit_report.md"


class StrictAlphaExistenceBitVsUnitOrientationSelectorAuditProbeTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        subprocess.run([sys.executable, str(SCRIPT)], check=True, cwd=ROOT)
        cls.payload = json.loads(REPORT_JSON.read_text(encoding="utf-8"))

    def test_result_kind_and_model(self):
        payload = self.payload
        self.assertEqual(
            payload["result_kind"],
            "SCRATCH_STRICT_ALPHA_EXISTENCE_BIT_VS_UNIT_ORIENTATION_SELECTOR_AUDIT_PROBE__NOT_A_THEOREM",
        )
        self.assertEqual(payload["status"], "existence-bit-can-gate-unoriented-orbit-not-select-d5-unit-branch")
        model = payload["finite_model"]
        self.assertEqual(model["ring"], "Z_12")
        self.assertEqual(model["support_count"], 792)
        self.assertEqual(model["histogram_class_count"], 35)
        self.assertEqual(model["automorphism_units"], [1, 5, 7, 11])
        self.assertEqual(model["existence_states"], [0, 1])
        self.assertEqual(model["target_q_power"], "256/243")
        self.assertEqual(model["target_eta"], "9/5")

    def test_bit_character_table(self):
        rows = {row["datum"]: row for row in self.payload["bit_character_table"]}
        existence = rows["existence_bit_E"]
        self.assertEqual(existence["aut_action"], "trivial/fixed")
        self.assertTrue(existence["can_gate_nonempty_state"])
        self.assertFalse(existence["can_select_A5_over_A1_without_extra_bit"])
        self.assertEqual(existence["symmetry_type"], "Aut-invariant scalar")

        unit = rows["unit_orientation_bit_U"]
        self.assertIn("unit mirror", unit["aut_action"])
        self.assertTrue(unit["can_select_A5_over_A1_without_extra_bit"])
        self.assertEqual(unit["symmetry_type"], "Aut-breaking character")

    def test_existence_gate_function_audit(self):
        audit = self.payload["existence_gate_function_audit"]
        self.assertEqual(audit["function_count"], 64)
        self.assertEqual(audit["canonical_existence_gate_count"], 1)
        self.assertEqual(audit["full_aut_invariant_singleton_d5_gate_count"], 0)
        example = audit["canonical_existence_gate_example"]
        self.assertEqual(example["E0_output"], "nothingness_only")
        self.assertEqual(example["E1_output"], "existence_unoriented_A1_A5_orbit")
        self.assertTrue(example["full_aut_invariant"])
        self.assertTrue(example["honest_existence_gate"])
        self.assertFalse(example["selects_singleton_d5_when_E1"])
        forbidden = [row for row in audit["rows"] if row["E1_output"] == "singleton_A5_d5"]
        self.assertTrue(forbidden)
        self.assertTrue(all(not row["full_aut_invariant"] for row in forbidden))

    def test_proof_interpretation_and_guardrails(self):
        proof = self.payload["exact_proof_certificate"]
        self.assertIn("fixes E=0 and E=1", proof["existence_bit_action"])
        self.assertIn("swap A1", proof["unit_branch_action"])
        self.assertIn("cannot choose a member", proof["type_mismatch"])
        self.assertIn("E=0 -> empty", proof["allowed_equation"])
        self.assertIn("requires a second", proof["blocked_equation"])
        self.assertIn("right symmetry type", proof["information_reading"])

        interpretation = self.payload["interpretation"]
        self.assertIn("not the same symmetry type", interpretation["direct_correction"])
        self.assertIn("E=0 -> empty", interpretation["what_it_can_do"])
        self.assertIn("silently adding", interpretation["what_it_cannot_do"])
        self.assertIn("Being information is not enough", interpretation["information_principle"])

        ontology = self.payload["ontology_guardrail"]
        self.assertIn("nadsoliton itself", ontology["allowed_reading"])
        self.assertIn("No separate informational layer", ontology["forbidden_reading"])

        hard_limits = "\n".join(self.payload["hard_limits"])
        self.assertIn("No identity K_legacy_ont == K_strict_gate", hard_limits)
        self.assertIn("No theorem derives the existence bit", hard_limits)
        self.assertIn("does not discharge QW-2191", hard_limits)
        self.assertIn("No QW-2191 discharge", hard_limits)
        self.assertIn("No ToE closure", hard_limits)
        self.assertTrue(REPORT_MD.exists())


if __name__ == "__main__":
    unittest.main()
