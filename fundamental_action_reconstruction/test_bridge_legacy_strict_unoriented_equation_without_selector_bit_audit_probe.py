import json
import subprocess
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "fundamental_action_reconstruction" / "scratch" / "bridge_legacy_strict_unoriented_equation_without_selector_bit_audit_probe.py"
REPORT_JSON = ROOT / "fundamental_action_reconstruction" / "scratch" / "bridge_legacy_strict_unoriented_equation_without_selector_bit_audit_report.json"
REPORT_MD = ROOT / "fundamental_action_reconstruction" / "scratch" / "bridge_legacy_strict_unoriented_equation_without_selector_bit_audit_report.md"


class LegacyStrictUnorientedEquationWithoutSelectorBitAuditProbeTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        subprocess.run([sys.executable, str(SCRIPT)], check=True, cwd=ROOT)
        cls.payload = json.loads(REPORT_JSON.read_text(encoding="utf-8"))

    def test_result_kind_model_and_kernel_split(self):
        payload = self.payload
        self.assertEqual(
            payload["result_kind"],
            "SCRATCH_LEGACY_STRICT_UNORIENTED_EQUATION_WITHOUT_SELECTOR_BIT_AUDIT_PROBE__NOT_A_THEOREM",
        )
        self.assertEqual(
            payload["status"],
            "quotient-bridge-equation-possible-singleton-d5-bridge-blocked-without-unit-bit",
        )
        model = payload["finite_model"]
        self.assertEqual(model["ring"], "Z_12")
        self.assertEqual(model["support_count"], 792)
        self.assertEqual(model["histogram_class_count"], 35)
        self.assertEqual(model["automorphism_units"], [1, 5, 7, 11])
        self.assertEqual(model["target_q_power"], "256/243")
        self.assertEqual(model["target_eta"], "9/5")

        kernels = payload["kernel_forms_kept_separate"]
        self.assertIn("K_legacy_ont", kernels["legacy_kernel_form"])
        self.assertIn("K_strict_gate", kernels["strict_kernel_form"])
        self.assertFalse(kernels["identity_asserted"])

    def test_aut_invariant_output_subset_audit(self):
        audit = self.payload["aut_invariant_output_subset_audit"]
        self.assertEqual(audit["selector_subset_count"], 4)
        self.assertEqual(audit["full_aut_invariant_output_names"], ["empty_relation", "unoriented_orbit_A1_A5"])
        self.assertFalse(audit["singleton_d5_full_aut_invariant"])
        self.assertTrue(audit["unoriented_orbit_output_available"])
        by_name = {row["name"]: row for row in audit["rows"]}
        self.assertFalse(by_name["singleton_A5_d5"]["full_aut_invariant"])
        self.assertTrue(by_name["unoriented_orbit_A1_A5"]["can_be_output_of_unoriented_bridge_relation"])

    def test_bridge_equation_classification_and_conditional_bit(self):
        by_type = {row["equation_type"]: row for row in self.payload["bridge_equation_classification"]}
        quotient = by_type["quotient_or_relation_bridge_without_selector_bit"]
        self.assertTrue(quotient["well_formed_without_selector_bit"])
        self.assertFalse(quotient["selects_singleton_d5"])
        d5 = by_type["deterministic_strict_d5_branch_bridge_without_selector_bit"]
        self.assertFalse(d5["well_formed_without_selector_bit"])
        self.assertTrue(d5["selects_singleton_d5"])
        conditional = by_type["conditional_oriented_bridge_with_selector_bit"]
        self.assertIn("only conditionally", conditional["selects_singleton_d5"])

        bit = self.payload["conditional_oriented_equation_after_bit"]
        self.assertTrue(bit["requires_external_or_derived_bit"])
        self.assertIn("{1,11}", bit["bit_meaning"])
        branches = {row["selector_bit_b"]: row["oriented_branch"] for row in bit["branches"]}
        self.assertEqual(branches, {0: "A1_k1_contiguous", 1: "A5_k5_d5"})

    def test_missing_slots_proof_interpretation_and_guardrails(self):
        slots = {row["slot"]: row for row in self.payload["missing_bridge_slots"]}
        self.assertIn("alpha_geo", slots["amplitude_normalization"]["legacy_side"])
        self.assertIn("beta_tors", slots["damping_renormalization"]["legacy_side"])
        self.assertIn("omega", slots["phase_frequency_translation"]["legacy_side"])
        self.assertIn("cannot be eliminated", slots["unit_orientation_selector"]["without_selector_bit_status"])

        proof = self.payload["exact_proof_certificate"]
        self.assertIn("no singleton fixed point", proof["fixed_point_fact"])
        self.assertIn("whole orbit", proof["relation_map_permission"])
        self.assertIn("b=1", proof["selector_bit_role"])

        interpretation = self.payload["interpretation"]
        self.assertIn("unoriented/quotient bridge equation", interpretation["direct_answer"])
        self.assertIn("silently imports the missing bit", interpretation["what_cannot_be_done_now"])
        self.assertIn("does not remove the need", interpretation["matter_shannon_relation"])

        ontology = self.payload["ontology_guardrail"]
        self.assertIn("No separate informational layer", ontology["forbidden_reading"])

        hard_limits = "\n".join(self.payload["hard_limits"])
        self.assertIn("No identity K_legacy_ont == K_strict_gate", hard_limits)
        self.assertIn("No theorem derives amplitude normalization", hard_limits)
        self.assertIn("singleton d5 branch selection is blocked", hard_limits)
        self.assertIn("No QW-2191 discharge", hard_limits)
        self.assertIn("No ToE closure", hard_limits)
        self.assertTrue(REPORT_MD.exists())


if __name__ == "__main__":
    unittest.main()
