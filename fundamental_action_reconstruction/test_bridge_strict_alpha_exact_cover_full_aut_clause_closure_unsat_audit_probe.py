import json
import subprocess
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "fundamental_action_reconstruction" / "scratch" / "bridge_strict_alpha_exact_cover_full_aut_clause_closure_unsat_audit_probe.py"
REPORT_JSON = ROOT / "fundamental_action_reconstruction" / "scratch" / "bridge_strict_alpha_exact_cover_full_aut_clause_closure_unsat_audit_report.json"
REPORT_MD = ROOT / "fundamental_action_reconstruction" / "scratch" / "bridge_strict_alpha_exact_cover_full_aut_clause_closure_unsat_audit_report.md"


class StrictAlphaExactCoverFullAutClauseClosureUnsatAuditProbeTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        subprocess.run([sys.executable, str(SCRIPT)], check=True, cwd=ROOT)
        cls.payload = json.loads(REPORT_JSON.read_text(encoding="utf-8"))

    def test_result_kind_and_model(self):
        payload = self.payload
        self.assertEqual(
            payload["result_kind"],
            "SCRATCH_STRICT_ALPHA_EXACT_COVER_FULL_AUT_CLAUSE_CLOSURE_UNSAT_AUDIT_PROBE__NO_GO_NOT_A_THEOREM",
        )
        self.assertEqual(
            payload["status"],
            "full-aut-closure-of-d1-d6-exact-cover-is-unsat-chi_11-kernel-certificate-only",
        )
        model = payload["finite_model"]
        self.assertEqual(model["ring"], "Z_12")
        self.assertEqual(model["active_count"], 5)
        self.assertEqual(model["automorphism_units"], [1, 5, 7, 11])
        self.assertEqual(model["base_forbidden_shells"], [1, 6])
        self.assertEqual(model["target_q_power"], "256/243")
        self.assertEqual(model["target_eta"], "9/5")

    def test_unit_shell_action_and_closure(self):
        audit = self.payload["closure_audit"]
        by_unit = {row["unit"]: row for row in audit["unit_shell_action"]}
        self.assertEqual(by_unit[1]["image_of_d1"], 1)
        self.assertEqual(by_unit[11]["image_of_d1"], 1)
        self.assertEqual(by_unit[5]["image_of_d1"], 5)
        self.assertEqual(by_unit[7]["image_of_d1"], 5)
        self.assertEqual({row["image_of_d6"] for row in by_unit.values()}, {6})
        self.assertTrue(by_unit[1]["preserves_base_clause_set_d1_d6"])
        self.assertTrue(by_unit[11]["preserves_base_clause_set_d1_d6"])
        self.assertFalse(by_unit[5]["preserves_base_clause_set_d1_d6"])
        self.assertFalse(by_unit[7]["preserves_base_clause_set_d1_d6"])
        self.assertEqual(audit["shell_orbit_of_d1_under_full_Aut"], [1, 5])
        self.assertEqual(audit["shell_orbit_of_d6_under_full_Aut"], [6])
        self.assertEqual(audit["full_Aut_closed_forbidden_shells"], [1, 5, 6])
        self.assertEqual(audit["chi_11_kernel_closed_forbidden_shells"], [1, 6])

    def test_base_exact_cover_survives_only_chi_11_kernel(self):
        base = self.payload["closure_audit"]["base_chi_11_system"]
        self.assertEqual(base["name"], "base_or_chi_11_closed_d1_plus_d6")
        self.assertEqual(base["forbidden_shells"], [1, 6])
        self.assertEqual(base["cardinality_compatible_assignments_checked"], 792)
        self.assertEqual(base["clauses_forbidden_pair_count"], 18)
        self.assertEqual(base["solution_count"], 12)
        self.assertEqual(base["dihedral_orbit_count"], 1)
        self.assertEqual(base["distance_histograms_d1_to_d6"], [[0, 3, 2, 1, 4, 0]])
        self.assertEqual(base["gap_necklaces"], [[2, 2, 3, 2, 3]])
        self.assertTrue(base["selects_A5_d5_histogram"])
        self.assertFalse(base["is_unsat"])

    def test_full_aut_closure_is_unsat(self):
        full = self.payload["closure_audit"]["full_Aut_closed_system"]
        self.assertEqual(full["name"], "full_Aut_closed_d1_plus_d5_plus_d6")
        self.assertEqual(full["forbidden_shells"], [1, 5, 6])
        self.assertEqual(full["cardinality_compatible_assignments_checked"], 792)
        self.assertEqual(full["clauses_forbidden_pair_count"], 30)
        self.assertEqual(full["solution_count"], 0)
        self.assertEqual(full["dihedral_orbit_count"], 0)
        self.assertEqual(full["distance_histograms_d1_to_d6"], [])
        self.assertEqual(full["gap_necklaces"], [])
        self.assertFalse(full["selects_A5_d5_histogram"])
        self.assertTrue(full["is_unsat"])

    def test_aut_conjugate_d5_plus_d6_selects_contiguous_side(self):
        conjugate = self.payload["closure_audit"]["aut_conjugate_d5_plus_d6_system"]
        self.assertEqual(conjugate["name"], "aut_conjugate_d5_plus_d6")
        self.assertEqual(conjugate["forbidden_shells"], [5, 6])
        self.assertEqual(conjugate["solution_count"], 12)
        self.assertEqual(conjugate["dihedral_orbit_count"], 1)
        self.assertEqual(conjugate["distance_histograms_d1_to_d6"], [[4, 3, 2, 1, 0, 0]])
        self.assertEqual(conjugate["gap_necklaces"], [[1, 1, 1, 1, 8]])
        self.assertFalse(conjugate["selects_A5_d5_histogram"])
        self.assertFalse(conjugate["is_unsat"])

    def test_proof_interpretation_and_guardrails(self):
        proof = self.payload["exact_proof_certificate"]
        self.assertIn("C(12,5)=792", proof["finite_domain"])
        self.assertIn("12 solutions", proof["base_certificate"])
        self.assertIn("orbit of d1 is {d1,d5}", proof["unit_closure_step"])
        self.assertIn("0 solutions", proof["unsat_result"])
        self.assertIn("chi_11-kernel", proof["conditional_positive"])
        self.assertIn("internal source", proof["source_obstruction"])

        interpretation = self.payload["interpretation"]
        self.assertIn("computationally correct", interpretation["what_was_proved"])
        self.assertIn("UNSAT", interpretation["new_no_go"])
        self.assertIn("Boolean clauses", interpretation["relation_to_previous_probe"])
        self.assertIn("QW-2191 remains open", interpretation["honest_limit"])

        ontology = self.payload["ontology_guardrail"]
        self.assertIn("nadsoliton itself", ontology["allowed_reading"])
        self.assertIn("No separate informational layer", ontology["forbidden_reading"])

        hard_limits = "\n".join(self.payload["hard_limits"])
        self.assertIn("No identity K_legacy_ont == K_strict_gate", hard_limits)
        self.assertIn("No theorem derives the chi_11-kernel", hard_limits)
        self.assertIn("finite full-Aut clause-closure UNSAT audit", hard_limits)
        self.assertIn("No QW-2191 discharge", hard_limits)
        self.assertIn("No ToE closure", hard_limits)
        self.assertTrue(REPORT_MD.exists())


if __name__ == "__main__":
    unittest.main()
