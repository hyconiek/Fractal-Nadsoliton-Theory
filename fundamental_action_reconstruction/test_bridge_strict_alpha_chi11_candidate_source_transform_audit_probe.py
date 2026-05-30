import json
import subprocess
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "fundamental_action_reconstruction" / "scratch" / "bridge_strict_alpha_chi11_candidate_source_transform_audit_probe.py"
REPORT_JSON = ROOT / "fundamental_action_reconstruction" / "scratch" / "bridge_strict_alpha_chi11_candidate_source_transform_audit_report.json"
REPORT_MD = ROOT / "fundamental_action_reconstruction" / "scratch" / "bridge_strict_alpha_chi11_candidate_source_transform_audit_report.md"


class StrictAlphaChi11CandidateSourceTransformAuditProbeTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        subprocess.run([sys.executable, str(SCRIPT)], check=True, cwd=ROOT)
        cls.payload = json.loads(REPORT_JSON.read_text(encoding="utf-8"))

    def test_result_kind_and_model(self):
        payload = self.payload
        self.assertEqual(
            payload["result_kind"],
            "SCRATCH_STRICT_ALPHA_CHI11_CANDIDATE_SOURCE_TRANSFORM_AUDIT_PROBE__NO_GO_NOT_A_THEOREM",
        )
        self.assertEqual(
            payload["status"],
            "candidate-records-do-not-export-chi11-without-shell-label-or-reduced-symmetry-import",
        )
        model = payload["finite_model"]
        self.assertEqual(model["ring"], "Z_12")
        self.assertEqual(model["active_count"], 5)
        self.assertEqual(model["automorphism_units"], [1, 5, 7, 11])
        self.assertEqual(model["branch_modes"], [1, 5, 7, 11])
        self.assertEqual(model["required_chi11_kernel"], [1, 11])
        self.assertEqual(model["required_chi11_nonzero_coset"], [5, 7])
        self.assertEqual(model["target_q_power"], "256/243")
        self.assertEqual(model["target_eta"], "9/5")

    def test_branch_and_unit_action_rows(self):
        branches = {row["name"]: row for row in self.payload["branch_rows"]}
        self.assertEqual(branches["A1_k1"]["branch_pair"], "contiguous_pair_A1_A11")
        self.assertEqual(branches["A11_k11"]["branch_pair"], "contiguous_pair_A1_A11")
        self.assertEqual(branches["A5_k5"]["branch_pair"], "d5_pair_A5_A7")
        self.assertEqual(branches["A7_k7"]["branch_pair"], "d5_pair_A5_A7")
        self.assertEqual(branches["A1_k1"]["required_chi11_value"], -1)
        self.assertEqual(branches["A5_k5"]["required_chi11_value"], 1)
        self.assertEqual(branches["A1_k1"]["full_aut_orbit_histogram_O15_O2_O3_O4_O6"], [4, 3, 2, 1, 0])
        self.assertEqual(branches["A5_k5"]["full_aut_orbit_histogram_O15_O2_O3_O4_O6"], [4, 3, 2, 1, 0])

        actions = {row["unit"]: row for row in self.payload["unit_action_rows"]}
        self.assertTrue(actions[1]["preserves_d5_pair"])
        self.assertTrue(actions[11]["preserves_d5_pair"])
        self.assertFalse(actions[5]["preserves_d5_pair"])
        self.assertFalse(actions[7]["preserves_d5_pair"])
        self.assertEqual(actions[5]["support_action"]["A1_k1"], "A5_k5")
        self.assertEqual(actions[7]["support_action"]["A5_k5"], "A11_k11")
        self.assertEqual(actions[5]["required_character_value"], -1)
        self.assertEqual(actions[11]["required_character_value"], 1)

    def test_candidate_summary(self):
        summary = self.payload["candidate_summary"]
        self.assertEqual(summary["candidate_count"], 7)
        self.assertEqual(summary["full_aut_invariant_candidate_count"], 3)
        self.assertEqual(summary["pair_distinguishing_candidate_count"], 4)
        self.assertEqual(summary["numeric_chi11_covariant_candidate_count"], 2)
        self.assertEqual(summary["allowed_strict_chi11_source_candidate_count"], 0)
        self.assertEqual(
            summary["chi11_covariant_but_imported_candidates"],
            ["shell_labelled_d5_minus_d1_count", "shell_labelled_energy_difference_(d5+d6)-(d1+d6)"],
        )

    def test_candidate_transform_rows(self):
        rows = {row["candidate"]: row for row in self.payload["candidate_transform_rows"]}
        self.assertTrue(rows["support_size"]["full_aut_invariant_on_four_branch_orbit"])
        self.assertFalse(rows["support_size"]["distinguishes_contiguous_pair_from_d5_pair"])
        self.assertTrue(rows["full_aut_orbit_histogram_O15_O2_O3_O4_O6"]["full_aut_invariant_on_four_branch_orbit"])
        self.assertFalse(rows["full_aut_orbit_histogram_O15_O2_O3_O4_O6"]["distinguishes_contiguous_pair_from_d5_pair"])

        gap = rows["dihedral_gap_necklace"]
        self.assertTrue(gap["distinguishes_contiguous_pair_from_d5_pair"])
        self.assertFalse(gap["full_aut_invariant_on_four_branch_orbit"])
        self.assertEqual(gap["requires_shell_label_or_reduced_symmetry"], "reduced cyclic-order premise")

        raw = rows["raw_distance_histogram_d1_to_d6"]
        self.assertTrue(raw["distinguishes_contiguous_pair_from_d5_pair"])
        self.assertFalse(raw["full_aut_invariant_on_four_branch_orbit"])
        self.assertTrue(raw["requires_shell_label_or_reduced_symmetry"])

        chi = rows["shell_labelled_d5_minus_d1_count"]
        self.assertEqual(chi["values_by_branch"], {"A11_k11": -4, "A1_k1": -4, "A5_k5": 4, "A7_k7": 4})
        self.assertTrue(chi["numeric_chi11_covariant_on_branches"])
        self.assertFalse(chi["exports_allowed_strict_chi11_source"])
        self.assertTrue(chi["requires_shell_label_or_reduced_symmetry"])

    def test_proof_interpretation_and_guardrails(self):
        proof = self.payload["exact_proof_certificate"]
        self.assertIn("four unit branches", proof["finite_domain"])
        self.assertIn("kernel {1,11}", proof["required_transform"])
        self.assertIn("[]", proof["invariant_no_go"])
        self.assertIn("shell_labelled_d5_minus_d1_count", proof["covariant_imports"])
        self.assertIn("Allowed strict chi_11 source candidates found without imported shell labels: 0", proof["allowed_source_count"])
        self.assertIn("does not prove", proof["boundary"])

        interpretation = self.payload["interpretation"]
        self.assertIn("right chi_11 covariance", interpretation["honest_positive"])
        self.assertIn("No audited full-Aut invariant", interpretation["honest_negative"])
        self.assertIn("source question", interpretation["relation_to_previous_probe"])

        ontology = self.payload["ontology_guardrail"]
        self.assertIn("nadsoliton itself", ontology["allowed_reading"])
        self.assertIn("No separate informational layer", ontology["forbidden_reading"])

        hard_limits = "\n".join(self.payload["hard_limits"])
        self.assertIn("No identity K_legacy_ont == K_strict_gate", hard_limits)
        self.assertIn("No theorem derives the chi_11-kernel", hard_limits)
        self.assertIn("candidate-source transform audit", hard_limits)
        self.assertIn("No QW-2191 discharge", hard_limits)
        self.assertIn("No ToE closure", hard_limits)
        self.assertTrue(REPORT_MD.exists())


if __name__ == "__main__":
    unittest.main()
