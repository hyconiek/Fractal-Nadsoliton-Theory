import json
import subprocess
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "fundamental_action_reconstruction" / "scratch" / "bridge_strict_alpha_full_aut_forbidden_shell_lattice_audit_probe.py"
REPORT_JSON = ROOT / "fundamental_action_reconstruction" / "scratch" / "bridge_strict_alpha_full_aut_forbidden_shell_lattice_audit_report.json"
REPORT_MD = ROOT / "fundamental_action_reconstruction" / "scratch" / "bridge_strict_alpha_full_aut_forbidden_shell_lattice_audit_report.md"


class StrictAlphaFullAutForbiddenShellLatticeAuditProbeTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        subprocess.run([sys.executable, str(SCRIPT)], check=True, cwd=ROOT)
        cls.payload = json.loads(REPORT_JSON.read_text(encoding="utf-8"))

    def test_result_kind_and_model(self):
        payload = self.payload
        self.assertEqual(
            payload["result_kind"],
            "SCRATCH_STRICT_ALPHA_FULL_AUT_FORBIDDEN_SHELL_LATTICE_AUDIT_PROBE__NO_GO_NOT_A_THEOREM",
        )
        self.assertEqual(
            payload["status"],
            "full-aut-forbidden-shell-lattice-classified-previous-d1-d5-d6-closure-is-minimal-but-not-unique-unsat-blocker",
        )
        model = payload["finite_model"]
        self.assertEqual(model["ring"], "Z_12")
        self.assertEqual(model["automorphism_units"], [1, 5, 7, 11])
        self.assertEqual(model["folded_shells"], [1, 2, 3, 4, 5, 6])
        self.assertEqual(model["computed_shell_orbits"], [[1, 5], [2], [3], [4], [6]])
        self.assertEqual(model["full_Aut_shell_orbits_used"], [[1, 5], [2], [3], [4], [6]])
        self.assertEqual(model["target_active_count"], 5)
        self.assertEqual(model["previous_full_Aut_closure"], [1, 5, 6])
        self.assertEqual(model["target_q_power"], "256/243")
        self.assertEqual(model["target_eta"], "9/5")

    def test_lattice_summary(self):
        summary = self.payload["lattice_summary"]
        self.assertEqual(summary["total_full_Aut_masks"], 32)
        self.assertEqual(summary["five_active_unsat_mask_count"], 25)
        self.assertEqual(summary["five_active_sat_mask_count"], 7)
        self.assertEqual(summary["minimal_five_unsat_mask_count"], 5)
        self.assertEqual(
            summary["minimal_five_unsat_forbidden_shells"],
            [[4], [2, 3], [3, 6], [1, 2, 5], [1, 5, 6]],
        )
        self.assertTrue(summary["previous_closure_is_minimal_five_unsat_blocker"])
        self.assertEqual(summary["previous_closure_alpha"], 3)
        self.assertEqual(summary["previous_closure_five_support_count"], 0)

    def test_minimal_unsat_rows(self):
        rows = self.payload["minimal_five_unsat_rows"]
        by_shells = {tuple(row["forbidden_shells"]): row for row in rows}
        self.assertEqual(set(by_shells), {(4,), (2, 3), (3, 6), (1, 2, 5), (1, 5, 6)})
        self.assertEqual(by_shells[(1, 5, 6)]["orbit_labels"], ["O15_d1_d5", "O6_d6"])
        self.assertEqual(by_shells[(1, 5, 6)]["maximum_independent_size"], 3)
        self.assertEqual(by_shells[(1, 5, 6)]["independence_profile_k0_to_k6"], [1, 12, 36, 16, 0, 0, 0])
        self.assertTrue(by_shells[(1, 5, 6)]["is_previous_full_Aut_closure"])
        self.assertEqual(by_shells[(4,)]["maximum_independent_size"], 4)
        self.assertEqual(by_shells[(2, 3)]["five_support_count"], 0)
        for row in rows:
            self.assertTrue(row["five_active_unsat"])
            self.assertEqual(row["five_support_count"], 0)

    def test_satisfiable_boundary_rows(self):
        sat_shells = [row["forbidden_shells"] for row in self.payload["five_active_satisfiable_rows"]]
        self.assertEqual(sat_shells, [[], [1, 5], [2], [3], [1, 3, 5], [6], [2, 6]])
        by_shells = {tuple(row["forbidden_shells"]): row for row in self.payload["five_active_satisfiable_rows"]}
        self.assertEqual(by_shells[(1, 5)]["five_support_count"], 12)
        self.assertEqual(by_shells[(6,)]["five_support_count"], 192)
        self.assertEqual(by_shells[(2, 6)]["five_support_count"], 24)
        for row in self.payload["five_active_satisfiable_rows"]:
            self.assertFalse(row["five_active_unsat"])
            self.assertGreater(row["five_support_count"], 0)

    def test_all_lattice_rows_and_previous_closure(self):
        rows = self.payload["all_lattice_rows"]
        self.assertEqual(len(rows), 32)
        previous_rows = [row for row in rows if row["is_previous_full_Aut_closure"]]
        self.assertEqual(len(previous_rows), 1)
        previous = previous_rows[0]
        self.assertEqual(previous["mask"], 17)
        self.assertEqual(previous["forbidden_shells"], [1, 5, 6])
        self.assertEqual(previous["maximum_independent_size"], 3)
        self.assertEqual(previous["five_support_count"], 0)
        self.assertEqual(previous["five_gap_necklaces"], [])

    def test_proof_interpretation_and_guardrails(self):
        proof = self.payload["exact_proof_certificate"]
        self.assertIn("exactly 32", proof["full_aut_lattice"])
        self.assertIn("25 make k=5 UNSAT", proof["classification_counts"])
        self.assertIn("5 inclusion-minimal", proof["minimal_blockers"])
        self.assertIn("[1, 5, 6] is one minimal blocker", proof["previous_closure_position"])
        self.assertIn("does not derive", proof["non_uniqueness_warning"])
        self.assertIn("d1-vs-d5 shell-label premise", proof["conditional_selector_reading"])

        interpretation = self.payload["interpretation"]
        self.assertIn("exclusion lattice", interpretation["what_was_added"])
        self.assertIn("inclusion-minimal", interpretation["honest_positive"])
        self.assertIn("not the unique", interpretation["honest_negative"])

        ontology = self.payload["ontology_guardrail"]
        self.assertIn("nadsoliton itself", ontology["allowed_reading"])
        self.assertIn("No separate informational layer", ontology["forbidden_reading"])

        hard_limits = "\n".join(self.payload["hard_limits"])
        self.assertIn("No identity K_legacy_ont == K_strict_gate", hard_limits)
        self.assertIn("No theorem derives the chi_11-kernel", hard_limits)
        self.assertIn("finite full-Aut forbidden-shell lattice audit", hard_limits)
        self.assertIn("No QW-2191 discharge", hard_limits)
        self.assertIn("No ToE closure", hard_limits)
        self.assertTrue(REPORT_MD.exists())


if __name__ == "__main__":
    unittest.main()
