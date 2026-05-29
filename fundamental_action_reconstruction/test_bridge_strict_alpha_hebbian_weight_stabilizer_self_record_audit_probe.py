import json
import subprocess
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "fundamental_action_reconstruction" / "scratch" / "bridge_strict_alpha_hebbian_weight_stabilizer_self_record_audit_probe.py"
REPORT_JSON = ROOT / "fundamental_action_reconstruction" / "scratch" / "bridge_strict_alpha_hebbian_weight_stabilizer_self_record_audit_report.json"
REPORT_MD = ROOT / "fundamental_action_reconstruction" / "scratch" / "bridge_strict_alpha_hebbian_weight_stabilizer_self_record_audit_report.md"


class StrictAlphaHebbianWeightStabilizerSelfRecordAuditProbeTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        subprocess.run([sys.executable, str(SCRIPT)], check=True, cwd=ROOT)
        cls.payload = json.loads(REPORT_JSON.read_text(encoding="utf-8"))

    def test_result_kind_and_model(self):
        payload = self.payload
        self.assertEqual(
            payload["result_kind"],
            "SCRATCH_STRICT_ALPHA_HEBBIAN_WEIGHT_STABILIZER_SELF_RECORD_AUDIT_PROBE__NOT_A_THEOREM",
        )
        self.assertEqual(
            payload["status"],
            "d5-hebbian-weight-carries-required-subgroup-conditionally-not-origin-theorem",
        )
        model = payload["finite_model"]
        self.assertEqual(model["ring"], "Z_12")
        self.assertEqual(model["active_count"], 5)
        self.assertEqual(model["automorphism_units"], [1, 5, 7, 11])
        self.assertEqual(model["minimal_required_subgroup_from_previous_probe"], [1, 11])
        self.assertEqual(model["target_q_power"], "256/243")
        self.assertEqual(model["target_eta"], "9/5")

    def test_d5_weight_self_record_certificate(self):
        certificate = self.payload["d5_weight_self_record_certificate"]
        self.assertEqual(certificate["unit_stabilizer"], [1, 11])
        self.assertTrue(certificate["equals_minimal_required_subgroup"])
        self.assertEqual(
            certificate["first_row_by_folded_distance_1_to_6"],
            {"1": "-25/12", "2": "11/12", "3": "-1/12", "4": "-13/12", "5": "23/12", "6": "-25/12"},
        )
        self.assertEqual(certificate["unit_5_obstruction_example"]["source_folded_distance"], 1)
        self.assertEqual(certificate["unit_5_obstruction_example"]["image_folded_distance"], 5)
        self.assertEqual(certificate["unit_5_obstruction_example"]["source_weight"], "-25/12")
        self.assertEqual(certificate["unit_5_obstruction_example"]["image_weight"], "23/12")
        self.assertIn("preserves exactly units {1,11}", certificate["readout"])

    def test_cross_teacher_nonuniqueness(self):
        audit = self.payload["cross_teacher_stabilizer_audit"]
        self.assertEqual(
            audit["classes_with_minimal_required_stabilizer"],
            ["contiguous_step_1_11", "fifth_step_d5_step_5_7"],
        )
        self.assertEqual(
            audit["classes_with_full_aut_stabilizer"],
            ["fourth_step_4_8_degenerate", "nyquist_step_6", "parity_minus_one_step_2_10", "third_step_3_9_degenerate"],
        )
        self.assertTrue(audit["contiguous_has_same_stabilizer_as_d5"])
        self.assertIn("not unique to d5", audit["nonuniqueness_warning"])
        d5 = audit["class_certificates"]["fifth_step_d5_step_5_7"]
        unit5 = next(row for row in d5["unit_obstruction_rows"] if row["unit"] == 5)
        self.assertFalse(unit5["preserves_weight_matrix"])
        self.assertEqual(unit5["mismatch_count_on_first_row_folded_distances"], 2)

    def test_interpretation_and_guardrails(self):
        payload = self.payload
        interpretation = payload["candidate_source_interpretation"]
        self.assertIn("same subgroup needed", interpretation["finite_gain"])
        self.assertIn("learned state itself", interpretation["conditional_positive_result"])
        self.assertIn("does not derive why the trace is d5", interpretation["honest_limit"])
        self.assertIn("cycle-metric probe", interpretation["relation_to_cycle_metric_probe"])

        ontology = payload["ontology_guardrail"]
        self.assertIn("nadsoliton itself", ontology["allowed_reading"])
        self.assertIn("No separate informational layer", ontology["forbidden_reading"])

        hard_limits = "\n".join(payload["hard_limits"])
        self.assertIn("No identity K_legacy_ont == K_strict_gate", hard_limits)
        self.assertIn("No theorem derives the d5 teacher/self-record trace", hard_limits)
        self.assertIn("No theorem derives a Hebbian learning law", hard_limits)
        self.assertIn("not a d5 origin theorem", hard_limits)
        self.assertIn("No QW-2191 discharge", hard_limits)
        self.assertIn("No ToE closure", hard_limits)
        self.assertTrue(REPORT_MD.exists())


if __name__ == "__main__":
    unittest.main()
