import json
import subprocess
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "fundamental_action_reconstruction" / "scratch" / "bridge_strict_alpha_full_aut_gap_independence_certificate_probe.py"
REPORT_JSON = ROOT / "fundamental_action_reconstruction" / "scratch" / "bridge_strict_alpha_full_aut_gap_independence_certificate_report.json"
REPORT_MD = ROOT / "fundamental_action_reconstruction" / "scratch" / "bridge_strict_alpha_full_aut_gap_independence_certificate_report.md"


class StrictAlphaFullAutGapIndependenceCertificateProbeTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        subprocess.run([sys.executable, str(SCRIPT)], check=True, cwd=ROOT)
        cls.payload = json.loads(REPORT_JSON.read_text(encoding="utf-8"))

    def test_result_kind_and_model(self):
        payload = self.payload
        self.assertEqual(
            payload["result_kind"],
            "SCRATCH_STRICT_ALPHA_FULL_AUT_GAP_INDEPENDENCE_CERTIFICATE_PROBE__NO_GO_NOT_A_THEOREM",
        )
        self.assertEqual(
            payload["status"],
            "full-aut-d1-d5-d6-closure-has-independence-number-3-so-five-active-exact-cover-is-unsat",
        )
        model = payload["finite_model"]
        self.assertEqual(model["ring"], "Z_12")
        self.assertEqual(model["automorphism_units"], [1, 5, 7, 11])
        self.assertEqual(model["forbidden_shells_full_Aut_closed"], [1, 5, 6])
        self.assertEqual(model["target_active_count"], 5)
        self.assertEqual(model["target_q_power"], "256/243")
        self.assertEqual(model["target_eta"], "9/5")

    def test_independence_profile_and_alpha(self):
        profile = {row["active_count"]: row for row in self.payload["independence_profile"]}
        self.assertEqual([profile[k]["independent_support_count"] for k in range(6)], [1, 12, 36, 16, 0, 0])
        self.assertTrue(profile[3]["exists"])
        self.assertFalse(profile[4]["exists"])
        self.assertFalse(profile[5]["exists"])
        self.assertEqual(profile[3]["gap_necklaces"], [[2, 2, 8], [4, 4, 4]])
        self.assertEqual(profile[4]["gap_necklaces"], [])
        self.assertEqual(profile[5]["gap_necklaces"], [])
        self.assertEqual(self.payload["maximum_independent_size"], 3)
        self.assertTrue(self.payload["target_five_active_is_unsat"])

    def test_gap_elimination_certificate_for_k4(self):
        certificate = self.payload["gap_elimination_certificate"]
        k4_rows = certificate["k4_rows"]
        self.assertEqual(len(k4_rows), 8)
        self.assertTrue(certificate["all_k4_necklaces_eliminated_by_d5_or_d6"])
        self.assertEqual(
            [row["gap_necklace"] for row in k4_rows],
            [
                [2, 2, 2, 6],
                [2, 2, 3, 5],
                [2, 2, 4, 4],
                [2, 3, 2, 5],
                [2, 3, 3, 4],
                [2, 3, 4, 3],
                [2, 4, 2, 4],
                [3, 3, 3, 3],
            ],
        )
        for row in k4_rows:
            self.assertTrue(row["all_gaps_at_least_2_from_d1_ban"])
            self.assertFalse(row["survives_d5_d6_bans"])
            self.assertIn(row["forbidden_interval_certificate"]["folded_distance"], [5, 6])

    def test_gap_elimination_certificate_for_k5(self):
        certificate = self.payload["gap_elimination_certificate"]
        k5_rows = certificate["k5_rows"]
        self.assertEqual(len(k5_rows), 3)
        self.assertTrue(certificate["all_k5_necklaces_eliminated_by_d5_or_d6"])
        self.assertEqual(
            [row["gap_necklace"] for row in k5_rows],
            [[2, 2, 2, 2, 4], [2, 2, 2, 3, 3], [2, 2, 3, 2, 3]],
        )
        for row in k5_rows:
            self.assertTrue(row["all_gaps_at_least_2_from_d1_ban"])
            self.assertFalse(row["survives_d5_d6_bans"])
            self.assertIn(row["forbidden_interval_certificate"]["folded_distance"], [5, 6])

    def test_proof_interpretation_and_guardrails(self):
        proof = self.payload["exact_proof_certificate"]
        self.assertIn("forbidden edges at folded distances d1,d5,d6", proof["graph_definition"])
        self.assertIn("d1 ban forces every cyclic gap", proof["gap_reduction"])
        self.assertIn("All 8", proof["four_support_elimination"])
        self.assertIn("All 3", proof["five_support_elimination"])
        self.assertIn("alpha(G)=3", proof["tightness"])
        self.assertIn("chi_11-conditional", proof["selector_consequence"])

        interpretation = self.payload["interpretation"]
        self.assertIn("alpha(G)=3", interpretation["what_was_added"])
        self.assertIn("8 k=4 and 3 k=5", interpretation["why_more_proof_like"])
        self.assertIn("QW-2191 discharge", interpretation["honest_limit"])

        ontology = self.payload["ontology_guardrail"]
        self.assertIn("nadsoliton itself", ontology["allowed_reading"])
        self.assertIn("No separate informational layer", ontology["forbidden_reading"])

        hard_limits = "\n".join(self.payload["hard_limits"])
        self.assertIn("No identity K_legacy_ont == K_strict_gate", hard_limits)
        self.assertIn("No theorem derives the chi_11-kernel", hard_limits)
        self.assertIn("finite full-Aut graph/gap no-go certificate", hard_limits)
        self.assertIn("No QW-2191 discharge", hard_limits)
        self.assertIn("No ToE closure", hard_limits)
        self.assertTrue(REPORT_MD.exists())


if __name__ == "__main__":
    unittest.main()
