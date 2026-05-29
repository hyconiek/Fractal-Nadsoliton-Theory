import json
import subprocess
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "fundamental_action_reconstruction" / "scratch" / "bridge_strict_alpha_unit_orbit_pair_selector_no_go_audit_probe.py"
REPORT_JSON = ROOT / "fundamental_action_reconstruction" / "scratch" / "bridge_strict_alpha_unit_orbit_pair_selector_no_go_audit_report.json"
REPORT_MD = ROOT / "fundamental_action_reconstruction" / "scratch" / "bridge_strict_alpha_unit_orbit_pair_selector_no_go_audit_report.md"


class StrictAlphaUnitOrbitPairSelectorNoGoAuditProbeTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        subprocess.run([sys.executable, str(SCRIPT)], check=True, cwd=ROOT)
        cls.payload = json.loads(REPORT_JSON.read_text(encoding="utf-8"))

    def test_result_kind_and_model(self):
        payload = self.payload
        self.assertEqual(
            payload["result_kind"],
            "SCRATCH_STRICT_ALPHA_UNIT_ORBIT_PAIR_SELECTOR_NO_GO_AUDIT_PROBE__NO_GO_NOT_A_THEOREM",
        )
        self.assertEqual(
            payload["status"],
            "full-aut-on-branch-pairs-forbids-singleton-d5-selector-chi_11-kernel-allows-it-conditionally",
        )
        model = payload["finite_model"]
        self.assertEqual(model["ring"], "Z_12")
        self.assertEqual(model["active_count"], 5)
        self.assertEqual(model["automorphism_units"], [1, 5, 7, 11])
        self.assertEqual(model["branch_modes"], [1, 5, 7, 11])
        self.assertEqual(model["target_q_power"], "256/243")
        self.assertEqual(model["target_eta"], "9/5")

    def test_branch_supports_and_pairs(self):
        by_name = {row["name"]: row for row in self.payload["branch_supports"]}
        self.assertEqual(by_name["A1_k1"]["support"], [0, 1, 2, 3, 4])
        self.assertEqual(by_name["A1_k1"]["distance_histogram_d1_to_d6"], [4, 3, 2, 1, 0, 0])
        self.assertEqual(by_name["A1_k1"]["orientation_pair"], "contiguous_pair_A1_A11")
        self.assertEqual(by_name["A11_k11"]["orientation_pair"], "contiguous_pair_A1_A11")
        self.assertEqual(by_name["A5_k5"]["support"], [0, 3, 5, 8, 10])
        self.assertEqual(by_name["A5_k5"]["distance_histogram_d1_to_d6"], [0, 3, 2, 1, 4, 0])
        self.assertEqual(by_name["A5_k5"]["orientation_pair"], "d5_pair_A5_A7")
        self.assertEqual(by_name["A7_k7"]["orientation_pair"], "d5_pair_A5_A7")

    def test_unit_action_on_pairs(self):
        by_unit = {row["unit"]: row for row in self.payload["unit_action_on_supports_and_pairs"]}
        self.assertEqual(by_unit[1]["pair_action"]["d5_pair_A5_A7"], "d5_pair_A5_A7")
        self.assertEqual(by_unit[11]["pair_action"]["d5_pair_A5_A7"], "d5_pair_A5_A7")
        self.assertEqual(by_unit[5]["pair_action"]["d5_pair_A5_A7"], "contiguous_pair_A1_A11")
        self.assertEqual(by_unit[7]["pair_action"]["contiguous_pair_A1_A11"], "d5_pair_A5_A7")
        self.assertTrue(by_unit[1]["preserves_d5_pair"])
        self.assertTrue(by_unit[11]["preserves_d5_pair"])
        self.assertFalse(by_unit[5]["preserves_d5_pair"])
        self.assertFalse(by_unit[7]["preserves_d5_pair"])
        self.assertEqual(by_unit[5]["support_action"]["A1_k1"], "A5_k5")
        self.assertEqual(by_unit[11]["support_action"]["A5_k5"], "A7_k7")

    def test_selector_subset_audit(self):
        audit = self.payload["selector_subset_audit"]
        self.assertEqual(audit["pair_universe"], ["contiguous_pair_A1_A11", "d5_pair_A5_A7"])
        self.assertEqual(audit["full_Aut_singleton_d5_selector_count"], 0)
        self.assertEqual(audit["chi_11_kernel_singleton_d5_selector_count"], 1)
        self.assertEqual(
            [row["selector"] for row in audit["full_Aut_invariant_selectors"]],
            ["empty_selector", "both_pairs_unoriented_orbit"],
        )
        reduced_selectors = {row["selector"] for row in audit["chi_11_kernel_invariant_selectors"]}
        self.assertIn("singleton_d5_pair", reduced_selectors)
        singleton_d5_rows = [row for row in audit["all_selector_rows"] if row["selector"] == "singleton_d5_pair"]
        self.assertEqual(
            {row["invariant_under"]: row["is_invariant"] for row in singleton_d5_rows},
            {"full_Aut_Z12": False, "chi_11_kernel_units_{1,11}": True},
        )

    def test_proof_interpretation_and_guardrails(self):
        proof = self.payload["exact_proof_certificate"]
        self.assertIn("one full Aut", proof["unit_orbit"])
        self.assertIn("contiguous_pair", proof["orientation_pairs"])
        self.assertIn("Units 5 and 7 swap", proof["pair_action"])
        self.assertIn("singleton d5_pair is not invariant", proof["selector_no_go"])
        self.assertIn("{1,11}", proof["conditional_positive"])
        self.assertIn("missing selector premise", proof["source_obstruction"])

        interpretation = self.payload["interpretation"]
        self.assertIn("cannot select singleton d5_pair", interpretation["direct_no_go"])
        self.assertIn("chi_11 kernel", interpretation["conditional_positive"])
        self.assertIn("branch-pair action", interpretation["relation_to_previous_probe"])
        self.assertIn("QW-2191 remains open", interpretation["honest_limit"])

        ontology = self.payload["ontology_guardrail"]
        self.assertIn("nadsoliton itself", ontology["allowed_reading"])
        self.assertIn("No separate informational layer", ontology["forbidden_reading"])

        hard_limits = "\n".join(self.payload["hard_limits"])
        self.assertIn("No identity K_legacy_ont == K_strict_gate", hard_limits)
        self.assertIn("No theorem derives the chi_11-kernel", hard_limits)
        self.assertIn("finite branch-pair no-go", hard_limits)
        self.assertIn("No QW-2191 discharge", hard_limits)
        self.assertIn("No ToE closure", hard_limits)
        self.assertTrue(REPORT_MD.exists())


if __name__ == "__main__":
    unittest.main()
