import json
import subprocess
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "fundamental_action_reconstruction" / "scratch" / "bridge_strict_alpha_chi11_shell_exclusion_lattice_unique_a5_selector_probe.py"
REPORT_JSON = ROOT / "fundamental_action_reconstruction" / "scratch" / "bridge_strict_alpha_chi11_shell_exclusion_lattice_unique_a5_selector_report.json"
REPORT_MD = ROOT / "fundamental_action_reconstruction" / "scratch" / "bridge_strict_alpha_chi11_shell_exclusion_lattice_unique_a5_selector_report.md"


class StrictAlphaChi11ShellExclusionLatticeUniqueA5SelectorProbeTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        subprocess.run([sys.executable, str(SCRIPT)], check=True, cwd=ROOT)
        cls.payload = json.loads(REPORT_JSON.read_text(encoding="utf-8"))

    def test_result_kind_and_model(self):
        payload = self.payload
        self.assertEqual(
            payload["result_kind"],
            "SCRATCH_STRICT_ALPHA_CHI11_SHELL_EXCLUSION_LATTICE_UNIQUE_A5_SELECTOR_PROBE__CONDITIONAL_NOT_A_THEOREM",
        )
        self.assertEqual(
            payload["status"],
            "chi11-shell-exclusion-lattice-has-unique-a5-d5-selector-d1-d6-conditional-on-shell-label-premise",
        )
        model = payload["finite_model"]
        self.assertEqual(model["ring"], "Z_12")
        self.assertEqual(model["active_count"], 5)
        self.assertEqual(model["folded_shells"], [1, 2, 3, 4, 5, 6])
        self.assertEqual(model["chi11_kernel_units"], [1, 11])
        self.assertEqual(model["full_Aut_units_for_contrast"], [1, 5, 7, 11])
        self.assertEqual(model["target_A5_histogram"], [0, 3, 2, 1, 4, 0])
        self.assertEqual(model["target_A1_histogram"], [4, 3, 2, 1, 0, 0])
        self.assertEqual(model["target_q_power"], "256/243")
        self.assertEqual(model["target_eta"], "9/5")

    def test_shell_action_audit(self):
        audit = self.payload["shell_action_audit"]
        chi_rows = {row["unit"]: row for row in audit["chi11_kernel_shell_action"]}
        self.assertTrue(chi_rows[1]["fixes_all_folded_shell_labels"])
        self.assertTrue(chi_rows[11]["fixes_all_folded_shell_labels"])
        self.assertEqual(chi_rows[11]["shell_images"]["d1"], "d1")
        self.assertEqual(chi_rows[11]["shell_images"]["d5"], "d5")

        full_rows = {row["unit"]: row for row in audit["full_Aut_shell_action_for_contrast"]}
        self.assertFalse(full_rows[5]["fixes_all_folded_shell_labels"])
        self.assertFalse(full_rows[7]["fixes_all_folded_shell_labels"])
        self.assertEqual(full_rows[5]["shell_images"]["d1"], "d5")
        self.assertEqual(full_rows[7]["shell_images"]["d5"], "d1")

    def test_lattice_summary(self):
        summary = self.payload["lattice_summary"]
        self.assertEqual(summary["total_chi11_shell_masks"], 64)
        self.assertEqual(summary["a5_d5_selecting_mask_count"], 1)
        self.assertEqual(summary["a5_d5_selecting_forbidden_shells"], [[1, 6]])
        self.assertEqual(summary["a5_d5_selecting_masks"], [33])
        self.assertEqual(summary["a1_contiguous_selecting_mask_count"], 1)
        self.assertEqual(summary["a1_contiguous_selecting_forbidden_shells"], [[5, 6]])
        self.assertEqual(summary["unsat_mask_count"], 51)
        self.assertEqual(summary["satisfiable_mask_count"], 13)
        self.assertEqual(
            summary["minimal_unsat_forbidden_shells"],
            [[4], [1, 2], [2, 3], [2, 5], [3, 6], [1, 5, 6]],
        )

    def test_unique_a5_selector_row(self):
        rows = self.payload["a5_d5_selector_rows"]
        self.assertEqual(len(rows), 1)
        row = rows[0]
        self.assertEqual(row["mask"], 33)
        self.assertEqual(row["forbidden_shells"], [1, 6])
        self.assertEqual(row["forbidden_shell_labels"], ["d1", "d6"])
        self.assertEqual(row["solution_count"], 12)
        self.assertEqual(row["dihedral_orbit_count"], 1)
        self.assertEqual(row["histograms_d1_to_d6"], [[0, 3, 2, 1, 4, 0]])
        self.assertEqual(row["gap_necklaces"], [[2, 2, 3, 2, 3]])
        self.assertTrue(row["selects_A5_d5"])
        self.assertFalse(row["selects_A1_contiguous"])
        self.assertFalse(row["is_unsat"])

    def test_conjugate_a1_boundary_row(self):
        rows = self.payload["a1_contiguous_selector_rows"]
        self.assertEqual(len(rows), 1)
        row = rows[0]
        self.assertEqual(row["mask"], 48)
        self.assertEqual(row["forbidden_shells"], [5, 6])
        self.assertEqual(row["forbidden_shell_labels"], ["d5", "d6"])
        self.assertEqual(row["solution_count"], 12)
        self.assertEqual(row["dihedral_orbit_count"], 1)
        self.assertEqual(row["histograms_d1_to_d6"], [[4, 3, 2, 1, 0, 0]])
        self.assertEqual(row["gap_necklaces"], [[1, 1, 1, 1, 8]])
        self.assertFalse(row["selects_A5_d5"])
        self.assertTrue(row["selects_A1_contiguous"])

    def test_all_lattice_rows_and_proof_guardrails(self):
        self.assertEqual(len(self.payload["all_lattice_rows"]), 64)
        proof = self.payload["exact_proof_certificate"]
        self.assertIn("64 chi_11-stable", proof["finite_domain"])
        self.assertIn("fix every folded shell label", proof["chi11_shell_stability"])
        self.assertIn("Exactly 1 mask selects the A5/d5", proof["unique_a5_selector"])
        self.assertIn("Exactly 1 mask selects the A1/contiguous", proof["conjugate_boundary"])
        self.assertIn("51 masks with no 5-support", proof["unsat_boundary"])
        self.assertIn("chi_11/shell-label premise", proof["conditional_scope"])

        interpretation = self.payload["interpretation"]
        self.assertIn("unique finite selector", interpretation["honest_positive"])
        self.assertIn("does not derive", interpretation["honest_negative"])
        self.assertIn("Full Aut destroys", interpretation["relation_to_full_Aut_audits"])

        ontology = self.payload["ontology_guardrail"]
        self.assertIn("nadsoliton itself", ontology["allowed_reading"])
        self.assertIn("No separate informational layer", ontology["forbidden_reading"])

        hard_limits = "\n".join(self.payload["hard_limits"])
        self.assertIn("No identity K_legacy_ont == K_strict_gate", hard_limits)
        self.assertIn("No theorem derives the chi_11-kernel", hard_limits)
        self.assertIn("conditional finite uniqueness certificate", hard_limits)
        self.assertIn("No QW-2191 discharge", hard_limits)
        self.assertIn("No ToE closure", hard_limits)
        self.assertTrue(REPORT_MD.exists())


if __name__ == "__main__":
    unittest.main()
