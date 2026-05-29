import json
import subprocess
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "fundamental_action_reconstruction" / "scratch" / "bridge_strict_alpha_full_aut_invariant_shell_energy_no_selector_audit_probe.py"
REPORT_JSON = ROOT / "fundamental_action_reconstruction" / "scratch" / "bridge_strict_alpha_full_aut_invariant_shell_energy_no_selector_audit_report.json"
REPORT_MD = ROOT / "fundamental_action_reconstruction" / "scratch" / "bridge_strict_alpha_full_aut_invariant_shell_energy_no_selector_audit_report.md"


class StrictAlphaFullAutInvariantShellEnergyNoSelectorAuditProbeTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        subprocess.run([sys.executable, str(SCRIPT)], check=True, cwd=ROOT)
        cls.payload = json.loads(REPORT_JSON.read_text(encoding="utf-8"))

    def test_result_kind_and_model(self):
        payload = self.payload
        self.assertEqual(
            payload["result_kind"],
            "SCRATCH_STRICT_ALPHA_FULL_AUT_INVARIANT_SHELL_ENERGY_NO_SELECTOR_AUDIT_PROBE__NO_GO_NOT_A_THEOREM",
        )
        self.assertEqual(payload["status"], "full-aut-invariant-linear-shell-energies-cannot-distinguish-A5-from-A1")
        model = payload["finite_model"]
        self.assertEqual(model["ring"], "Z_12")
        self.assertEqual(model["support_count"], 792)
        self.assertEqual(model["histogram_class_count"], 35)
        self.assertEqual(model["automorphism_units"], [1, 5, 7, 11])
        self.assertEqual(model["target_q_power"], "256/243")
        self.assertEqual(model["target_eta"], "9/5")
        self.assertEqual(model["A1_histogram_d1_to_d6"], [4, 3, 2, 1, 0, 0])
        self.assertEqual(model["A5_histogram_d1_to_d6"], [0, 3, 2, 1, 4, 0])

    def test_invariant_weight_structure(self):
        audit = self.payload["invariant_weight_structure_audit"]
        self.assertEqual(audit["shell_orbits_under_full_Aut"], [[1, 5], [2], [3], [4], [6]])
        self.assertEqual(audit["linear_invariance_condition"], "w1=w5; w2,w3,w4,w6 free")
        by_unit = {row["unit"]: row["shell_permutation_source_to_target"] for row in audit["unit_shell_action_rows"]}
        self.assertEqual(by_unit[1], [1, 2, 3, 4, 5, 6])
        self.assertEqual(by_unit[11], [1, 2, 3, 4, 5, 6])
        self.assertEqual(by_unit[5], [5, 2, 3, 4, 1, 6])
        self.assertEqual(by_unit[7], [5, 2, 3, 4, 1, 6])

        diff = audit["A1_A5_symbolic_score_difference"]
        self.assertEqual(diff["A1_minus_A5_histogram"], [4, 0, 0, 0, -4, 0])
        self.assertEqual(diff["dot_with_general_weight"], "4*w1 - 4*w5")
        self.assertEqual(diff["full_aut_invariance_condition"], "w1 = w5")
        self.assertTrue(diff["therefore_A1_A5_scores_equal"])

    def test_binary_mask_no_go(self):
        audit = self.payload["binary_full_aut_invariant_mask_minimizer_audit"]
        self.assertEqual(audit["binary_full_aut_invariant_mask_count"], 31)
        self.assertEqual(audit["unique_A5_selecting_mask_count"], 0)
        self.assertEqual(len(audit["rows"]), 31)
        self.assertTrue(all(not row["selects_unique_A5_d5_orbit"] for row in audit["rows"]))
        self.assertTrue(all(row["weight"][0] == row["weight"][4] for row in audit["rows"]))
        self.assertIn("d1_nearest + d4 + d5 + d6_antipodal", audit["mask_formulas_touching_A1_or_A5_winners"])

    def test_proof_interpretation_and_guardrails(self):
        proof = self.payload["exact_proof_certificate"]
        self.assertIn("swaps d1 with d5", proof["full_aut_shell_action"])
        self.assertIn("w1=w5", proof["invariant_weight_condition"])
        self.assertIn("zero A1/A5 score difference", proof["A1_A5_pair_obstruction"])
        self.assertIn("31 nonzero binary", proof["binary_enumeration"])
        self.assertIn("break full Aut", proof["source_obstruction"])

        interpretation = self.payload["interpretation"]
        self.assertIn("cannot distinguish", interpretation["direct_no_go"])
        self.assertIn("zero unique A5 selectors", interpretation["computational_confirmation"])
        self.assertIn("not full-Aut invariant", interpretation["relation_to_previous_probe"])
        self.assertIn("shell-linear", interpretation["honest_limit"])

        ontology = self.payload["ontology_guardrail"]
        self.assertIn("nadsoliton itself", ontology["allowed_reading"])
        self.assertIn("No separate informational layer", ontology["forbidden_reading"])

        hard_limits = "\n".join(self.payload["hard_limits"])
        self.assertIn("No identity K_legacy_ont == K_strict_gate", hard_limits)
        self.assertIn("No theorem derives a full-Aut invariant selector", hard_limits)
        self.assertIn("shell-linear no-go audit", hard_limits)
        self.assertIn("No QW-2191 discharge", hard_limits)
        self.assertIn("No ToE closure", hard_limits)
        self.assertTrue(REPORT_MD.exists())


if __name__ == "__main__":
    unittest.main()
