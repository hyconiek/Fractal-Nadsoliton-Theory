import json
import subprocess
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "fundamental_action_reconstruction" / "scratch" / "bridge_strict_alpha_hebbian_prospective_path_event_selector_no_go_probe.py"
REPORT_JSON = ROOT / "fundamental_action_reconstruction" / "scratch" / "bridge_strict_alpha_hebbian_prospective_path_event_selector_no_go_report.json"
REPORT_MD = ROOT / "fundamental_action_reconstruction" / "scratch" / "bridge_strict_alpha_hebbian_prospective_path_event_selector_no_go_report.md"


class StrictAlphaHebbianProspectivePathEventSelectorNoGoProbeTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        subprocess.run([sys.executable, str(SCRIPT)], check=True, cwd=ROOT)
        cls.payload = json.loads(REPORT_JSON.read_text(encoding="utf-8"))

    def test_result_kind_and_finite_model(self):
        payload = self.payload
        self.assertEqual(
            payload["result_kind"],
            "SCRATCH_STRICT_ALPHA_HEBBIAN_PROSPECTIVE_PATH_EVENT_SELECTOR_NO_GO_PROBE__NOT_A_THEOREM",
        )
        self.assertEqual(payload["status"], "aut-invariant-prospective-path-events-cannot-select-current-d5")
        model = payload["finite_model"]
        self.assertEqual(model["ring"], "Z_12")
        self.assertEqual(model["support_count"], 792)
        self.assertEqual(model["histogram_class_count"], 35)
        self.assertEqual(model["automorphism_units"], [1, 5, 7, 11])
        self.assertEqual(model["target_q_power"], "256/243")
        self.assertEqual(model["target_eta"], "9/5")
        self.assertEqual([axis["name"] for axis in model["survivor_axes"]], ["A1_k1_contiguous", "A5_k5_d5"])

    def test_path_horizon_audit_counts(self):
        audit = self.payload["path_event_horizon_audit"]
        self.assertEqual(audit["max_horizon"], 8)
        self.assertEqual(audit["enumerated_event_max_horizon"], 4)
        self.assertEqual(audit["enumerated_current_d5_dominant_total"], 0)
        self.assertEqual(audit["enumerated_terminal_d5_dominant_total"], 0)
        self.assertTrue(audit["all_exact_rows_have_noninvariant_current_A5_cylinder"])
        self.assertTrue(audit["all_exact_rows_have_noninvariant_terminal_A5_cylinder"])
        self.assertEqual([row["horizon"] for row in audit["rows"]], list(range(9)))

        for row in audit["rows"]:
            horizon = row["horizon"]
            self.assertEqual(row["path_length"], horizon + 1)
            self.assertEqual(row["path_count"], 2 ** (horizon + 1))
            self.assertEqual(row["mirror_orbit_count"], 2**horizon)
            self.assertEqual(row["fixed_path_count_under_mirror"], 0)
            self.assertEqual(row["global_unit_axis_bits_required"], 1)
            self.assertFalse(row["current_A5_cylinder_invariant"])
            self.assertFalse(row["terminal_A5_cylinder_invariant"])
            self.assertEqual(row["current_A5_cylinder_counts"]["current_A5_minus_A1"], 2**horizon)
            self.assertEqual(row["terminal_A5_cylinder_counts"]["terminal_A5_minus_A1"], 2**horizon)

    def test_enumerated_event_no_bias(self):
        rows = self.payload["path_event_horizon_audit"]["rows"][:5]
        expected_event_counts = [2, 4, 16, 256, 65536]
        for row, expected_count in zip(rows, expected_event_counts):
            self.assertEqual(row["enumerated_invariant_event_count"], expected_count)
            self.assertEqual(row["enumerated_current_d5_dominant_event_count"], 0)
            self.assertEqual(row["enumerated_terminal_d5_dominant_event_count"], 0)
            self.assertEqual(row["enumerated_nonzero_current_bias_count"], 0)
            self.assertEqual(row["enumerated_nonzero_terminal_bias_count"], 0)

    def test_sample_orbit_and_proof_certificate(self):
        sample = self.payload["sample_orbit_certificate"]
        self.assertEqual(sample["horizon"], 3)
        for orbit in sample["sample_orbits"]:
            self.assertTrue(orbit["orbit_has_one_current_A1_and_one_current_A5"])
            self.assertTrue(orbit["orbit_has_one_terminal_A1_and_one_terminal_A5"])
            self.assertEqual(len(orbit["paths"]), 2)

        proof = self.payload["exact_proof_certificate"]
        self.assertIn("Omega_H", proof["path_space"])
        self.assertIn("no fixed path", proof["mirror_action"])
        self.assertIn("union of whole mirror pairs", proof["invariant_event_condition"])
        self.assertIn("same number of paths", proof["current_balance_consequence"])
        self.assertIn("missing unit-axis bit", proof["selector_consequence"])

    def test_interpretation_and_guardrails(self):
        interpretation = self.payload["selector_interpretation"]
        self.assertIn("possible future paths", interpretation["question_tested"])
        self.assertIn("No full-Aut-invariant", interpretation["negative_result"])
        self.assertIn("one side of every mirror pair", interpretation["conditional_positive_result"])
        self.assertIn("highest-future-resonance", interpretation["relation_to_highest_future_resonance"])

        ontology = self.payload["ontology_guardrail"]
        self.assertIn("nadsoliton itself", ontology["allowed_reading"])
        self.assertIn("No separate informational layer", ontology["forbidden_reading"])

        hard_limits = "\n".join(self.payload["hard_limits"])
        self.assertIn("No identity K_legacy_ont == K_strict_gate", hard_limits)
        self.assertIn("No theorem derives a faith-like prospective path-record", hard_limits)
        self.assertIn("Full Aut(Z_12)-invariant finite future-path events", hard_limits)
        self.assertIn("No QW-2191 discharge", hard_limits)
        self.assertIn("No ToE closure", hard_limits)
        self.assertTrue(REPORT_MD.exists())


if __name__ == "__main__":
    unittest.main()
