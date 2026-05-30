import json
import subprocess
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "fundamental_action_reconstruction" / "scratch" / "bridge_strict_alpha_full_aut_block_amplitude_chi11_polarity_obstruction_probe.py"
REPORT_JSON = ROOT / "fundamental_action_reconstruction" / "scratch" / "bridge_strict_alpha_full_aut_block_amplitude_chi11_polarity_obstruction_report.json"
REPORT_MD = ROOT / "fundamental_action_reconstruction" / "scratch" / "bridge_strict_alpha_full_aut_block_amplitude_chi11_polarity_obstruction_report.md"


class StrictAlphaFullAutBlockAmplitudeChi11PolarityObstructionProbeTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        subprocess.run([sys.executable, str(SCRIPT)], check=True, cwd=ROOT)
        cls.payload = json.loads(REPORT_JSON.read_text(encoding="utf-8"))

    def test_result_kind_and_finite_model(self):
        payload = self.payload
        self.assertEqual(
            payload["result_kind"],
            "SCRATCH_STRICT_ALPHA_FULL_AUT_BLOCK_AMPLITUDE_CHI11_POLARITY_OBSTRUCTION__BLOCK_ONLY_NOT_A_SELECTOR_THEOREM",
        )
        self.assertEqual(
            payload["status"],
            "full-aut-amplitude-uniquely-locates-branch-block-but-not-chi11-polarity",
        )
        model = payload["finite_model"]
        self.assertEqual(model["ring"], "Z_12")
        self.assertEqual(model["active_count"], 5)
        self.assertEqual(model["support_count"], 792)
        self.assertEqual(model["dihedral_units"], [1, 11])
        self.assertEqual(model["automorphism_units"], [1, 5, 7, 11])
        self.assertEqual(model["d12_orbit_count"], 38)
        self.assertEqual(model["full_affine_aut_orbit_count"], 25)
        self.assertEqual(model["target_q_power"], "256/243")
        self.assertEqual(model["target_eta"], "9/5")

    def test_block_amplitude_summary_and_rows(self):
        summary = self.payload["block_amplitude_summary"]
        self.assertEqual(summary["candidate_score"], "max_component_abs(h_d5-h_d1) per full-Aut block")
        self.assertTrue(summary["full_aut_block_score"])
        self.assertFalse(summary["exports_chi11_polarity"])
        self.assertEqual(summary["max_amplitude"], 4)
        self.assertEqual(summary["maximizer_count"], 1)
        self.assertTrue(summary["unique_maximizer_is_branch_full_aut_block"])
        self.assertEqual(summary["amplitude_distribution"], {"0": 14, "1": 6, "2": 4, "4": 1})
        self.assertEqual(summary["branch_block_signed_values"], [-4, 4])
        self.assertTrue(summary["branch_block_has_polarity_pair"])

        rows = self.payload["full_aut_block_rows"]
        self.assertEqual(len(rows), 25)
        self.assertEqual(
            sorted(row["unoriented_abs_d5_minus_d1_amplitude"] for row in rows),
            [0] * 14 + [1] * 6 + [2] * 4 + [4],
        )
        self.assertEqual(sum(1 for row in rows if row["contains_branch_supports"]), 1)
        self.assertEqual(sum(1 for row in rows if row["unoriented_abs_d5_minus_d1_amplitude"] == 4), 1)

    def test_branch_block_and_branch_rows(self):
        block = self.payload["branch_block_certificate"]
        self.assertEqual(block["full_affine_orbit_index"], 0)
        self.assertEqual(block["dihedral_component_indices"], [0, 37])
        self.assertEqual(block["dihedral_component_count"], 2)
        self.assertEqual(block["signed_d5_minus_d1_values_by_component"], [-4, 4])
        self.assertEqual(block["unoriented_abs_d5_minus_d1_amplitude"], 4)
        self.assertTrue(block["has_polarity_pair"])
        self.assertTrue(block["contains_branch_supports"])
        self.assertEqual(block["component_rows"][0]["histogram_d1_to_d6"], [4, 3, 2, 1, 0, 0])
        self.assertEqual(block["component_rows"][1]["histogram_d1_to_d6"], [0, 3, 2, 1, 4, 0])

        branches = {row["name"]: row for row in self.payload["branch_rows"]}
        self.assertEqual(branches["A1_k1"]["full_affine_orbit_index"], 0)
        self.assertEqual(branches["A5_k5"]["full_affine_orbit_index"], 0)
        self.assertEqual(branches["A7_k7"]["full_affine_orbit_index"], 0)
        self.assertEqual(branches["A11_k11"]["full_affine_orbit_index"], 0)
        self.assertEqual(branches["A1_k1"]["d5_minus_d1"], -4)
        self.assertEqual(branches["A11_k11"]["d5_minus_d1"], -4)
        self.assertEqual(branches["A5_k5"]["d5_minus_d1"], 4)
        self.assertEqual(branches["A7_k7"]["d5_minus_d1"], 4)

    def test_proof_interpretation_and_guardrails(self):
        proof = self.payload["exact_proof_certificate"]
        self.assertIn("C(12,5)=792", proof["finite_domain"])
        self.assertIn("|h_d5-h_d1| is unchanged", proof["amplitude_invariance"])
        self.assertIn("max amplitude is 4", proof["unique_block_location"])
        self.assertIn("signed values -4 and +4", proof["polarity_obstruction"])
        self.assertIn("A5 over A1 still requires", proof["selector_boundary"])
        self.assertIn("does not discharge QW-2191", proof["strict_limit"])

        interpretation = self.payload["interpretation"]
        self.assertIn("identifies the full-Aut block", interpretation["honest_positive"])
        self.assertIn("polarity remains unavailable", interpretation["honest_negative"])
        self.assertIn("forgetting the sign", interpretation["relation_to_previous_probe"])

        ontology = self.payload["ontology_guardrail"]
        self.assertIn("nadsoliton itself", ontology["allowed_reading"])
        self.assertIn("No separate informational layer", ontology["forbidden_reading"])

        hard_limits = "\n".join(self.payload["hard_limits"])
        self.assertIn("No identity K_legacy_ont == K_strict_gate", hard_limits)
        self.assertIn("No theorem derives chi_11", hard_limits)
        self.assertIn("No QW-2191 discharge", hard_limits)
        self.assertIn("No ToE closure", hard_limits)
        self.assertIn("does not select chi_11 polarity", hard_limits)
        self.assertTrue(REPORT_MD.exists())


if __name__ == "__main__":
    unittest.main()
