import json
import subprocess
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "fundamental_action_reconstruction" / "scratch" / "bridge_strict_alpha_unit_character_cohomology_obstruction_audit_probe.py"
REPORT_JSON = ROOT / "fundamental_action_reconstruction" / "scratch" / "bridge_strict_alpha_unit_character_cohomology_obstruction_audit_report.json"
REPORT_MD = ROOT / "fundamental_action_reconstruction" / "scratch" / "bridge_strict_alpha_unit_character_cohomology_obstruction_audit_report.md"


class StrictAlphaUnitCharacterCohomologyObstructionAuditProbeTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        subprocess.run([sys.executable, str(SCRIPT)], check=True, cwd=ROOT)
        cls.payload = json.loads(REPORT_JSON.read_text(encoding="utf-8"))

    def test_result_kind_and_model(self):
        payload = self.payload
        self.assertEqual(
            payload["result_kind"],
            "SCRATCH_STRICT_ALPHA_UNIT_CHARACTER_COHOMOLOGY_OBSTRUCTION_AUDIT_PROBE__NOT_A_THEOREM",
        )
        self.assertEqual(
            payload["status"],
            "required-d5-unit-bit-is-nontrivial-character-not-invariant-scalar-coboundary",
        )
        model = payload["finite_model"]
        self.assertEqual(model["ring"], "Z_12")
        self.assertEqual(model["support_count"], 792)
        self.assertEqual(model["histogram_class_count"], 35)
        self.assertEqual(model["automorphism_units"], [1, 5, 7, 11])
        self.assertEqual(model["target_q_power"], "256/243")
        self.assertEqual(model["target_eta"], "9/5")

    def test_character_table_audit(self):
        audit = self.payload["character_table_audit"]
        self.assertEqual(audit["character_count"], 4)
        names = {row["name"] for row in audit["rows"]}
        self.assertEqual(names, {"trivial_character", "chi_11_required_d5_unit_axis", "chi_5_kernel", "chi_7_kernel"})
        required = audit["required_character"]
        self.assertEqual(required["name"], "chi_11_required_d5_unit_axis")
        self.assertEqual(required["kernel"], [1, 11])
        self.assertEqual(required["nonzero_coset"], [5, 7])
        self.assertTrue(required["is_required_d5_unit_axis_character"])

    def test_zero_cochain_coboundary_audit(self):
        audit = self.payload["zero_cochain_coboundary_audit"]
        self.assertEqual(audit["cochain_count"], 4)
        self.assertEqual(audit["invariant_scalar_cochain_count"], 2)
        self.assertEqual(audit["invariant_scalar_coboundary_required_chi_11_count"], 0)
        self.assertEqual(audit["noninvariant_cochain_required_chi_11_count"], 2)
        invariant_rows = [row for row in audit["rows"] if row["full_aut_invariant_scalar"]]
        self.assertTrue(all(row["coboundary_character_name"] == "trivial_character" for row in invariant_rows))
        noninvariant_rows = [row for row in audit["rows"] if not row["full_aut_invariant_scalar"]]
        self.assertTrue(all(row["coboundary_equals_required_chi_11"] for row in noninvariant_rows))
        self.assertTrue(all("imports the bit" in row["interpretation"] for row in noninvariant_rows))

    def test_subgroup_proof_interpretation_and_guardrails(self):
        subgroup = self.payload["subgroup_quotient_audit"]
        self.assertEqual(subgroup["required_subgroup"], [1, 11])
        self.assertEqual(subgroup["required_quotient_index"], 2)
        by_subgroup = {tuple(row["subgroup"]): row for row in subgroup["rows"]}
        self.assertTrue(by_subgroup[(1, 11)]["singleton_d5_invariant_under_subgroup"])
        self.assertTrue(by_subgroup[(1, 11)]["adjoins_required_character"])
        self.assertFalse(by_subgroup[(1, 5, 7, 11)]["singleton_d5_invariant_under_subgroup"])

        proof = self.payload["exact_proof_certificate"]
        self.assertIn("V4", proof["group"])
        self.assertIn("chi_11", proof["required_character"])
        self.assertIn("trivial character", proof["invariant_scalar_coboundary"])
        self.assertIn("already labelled", proof["noninvariant_cochain_warning"])
        self.assertIn("nontrivial Z2 character", proof["cohomology_reading"])

        interpretation = self.payload["interpretation"]
        self.assertIn("nontrivial unit character chi_11", interpretation["direct_result"])
        self.assertIn("Every full-Aut-invariant scalar", interpretation["negative_result"])
        self.assertIn("does not derive chi_11", interpretation["honest_limit"])

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
