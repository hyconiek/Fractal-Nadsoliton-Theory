import json
import subprocess
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "fundamental_action_reconstruction" / "scratch" / "bridge_strict_alpha_chi11_d1d6_local_necessity_certificate_probe.py"
REPORT_JSON = ROOT / "fundamental_action_reconstruction" / "scratch" / "bridge_strict_alpha_chi11_d1d6_local_necessity_certificate_report.json"
REPORT_MD = ROOT / "fundamental_action_reconstruction" / "scratch" / "bridge_strict_alpha_chi11_d1d6_local_necessity_certificate_report.md"


class StrictAlphaChi11D1D6LocalNecessityCertificateProbeTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        subprocess.run([sys.executable, str(SCRIPT)], check=True, cwd=ROOT)
        cls.payload = json.loads(REPORT_JSON.read_text(encoding="utf-8"))

    def test_result_kind_and_model(self):
        payload = self.payload
        self.assertEqual(
            payload["result_kind"],
            "SCRATCH_STRICT_ALPHA_CHI11_D1D6_LOCAL_NECESSITY_CERTIFICATE_PROBE__CONDITIONAL_NOT_A_THEOREM",
        )
        self.assertEqual(
            payload["status"],
            "d1-d6-is-locally-necessary-and-upward-isolated-inside-chi11-shell-exclusion-grammar",
        )
        model = payload["finite_model"]
        self.assertEqual(model["ring"], "Z_12")
        self.assertEqual(model["active_count"], 5)
        self.assertEqual(model["folded_shells"], [1, 2, 3, 4, 5, 6])
        self.assertEqual(model["chi11_kernel_units"], [1, 11])
        self.assertEqual(model["center_forbidden_shells"], [1, 6])
        self.assertEqual(model["target_A5_histogram"], [0, 3, 2, 1, 4, 0])
        self.assertEqual(model["target_A1_histogram"], [4, 3, 2, 1, 0, 0])
        self.assertEqual(model["target_q_power"], "256/243")
        self.assertEqual(model["target_eta"], "9/5")

    def test_local_boundary_summary(self):
        summary = self.payload["local_boundary_summary"]
        self.assertTrue(summary["center_selects_A5_d5_uniquely"])
        self.assertEqual(summary["deletion_neighbor_count"], 2)
        self.assertEqual(summary["addition_neighbor_count"], 4)
        self.assertTrue(summary["all_deletions_overselect"])
        self.assertTrue(summary["all_additions_unsat"])
        self.assertEqual(summary["delete_d1_solution_count"], 192)
        self.assertEqual(summary["delete_d6_solution_count"], 36)
        self.assertEqual(summary["addition_solution_counts"], {"add_d2": 0, "add_d3": 0, "add_d4": 0, "add_d5": 0})

    def test_center_and_deletion_neighbors(self):
        rows = {row["name"]: row for row in self.payload["local_boundary_rows"]}
        center = rows["center_d1_plus_d6"]
        self.assertEqual(center["forbidden_shells"], [1, 6])
        self.assertEqual(center["solution_count"], 12)
        self.assertEqual(center["dihedral_orbit_count"], 1)
        self.assertEqual(center["histogram_rows"], [{"histogram_d1_to_d6": [0, 3, 2, 1, 4, 0], "classification": "A5_d5_target"}])
        self.assertEqual(center["gap_necklaces"], [[2, 2, 3, 2, 3]])
        self.assertTrue(center["selects_A5_d5_uniquely"])

        delete_d1 = rows["delete_d1"]
        self.assertEqual(delete_d1["forbidden_shells"], [6])
        self.assertEqual(delete_d1["solution_count"], 192)
        self.assertEqual(delete_d1["histogram_class_count"], 9)
        self.assertEqual(delete_d1["gap_necklace_count"], 10)
        self.assertTrue(delete_d1["contains_A1_contiguous_histogram"])
        self.assertTrue(delete_d1["contains_A5_d5_histogram"])
        self.assertFalse(delete_d1["selects_A5_d5_uniquely"])

        delete_d6 = rows["delete_d6"]
        self.assertEqual(delete_d6["forbidden_shells"], [1])
        self.assertEqual(delete_d6["solution_count"], 36)
        self.assertEqual(delete_d6["histogram_class_count"], 3)
        self.assertEqual(delete_d6["gap_necklace_count"], 3)
        self.assertTrue(delete_d6["contains_A5_d5_histogram"])
        self.assertFalse(delete_d6["contains_A1_contiguous_histogram"])
        self.assertFalse(delete_d6["selects_A5_d5_uniquely"])

    def test_addition_neighbors_are_unsat(self):
        rows = {row["name"]: row for row in self.payload["local_boundary_rows"]}
        for name, forbidden_shells in {
            "add_d2": [1, 2, 6],
            "add_d3": [1, 3, 6],
            "add_d4": [1, 4, 6],
            "add_d5": [1, 5, 6],
        }.items():
            row = rows[name]
            self.assertEqual(row["forbidden_shells"], forbidden_shells)
            self.assertEqual(row["solution_count"], 0)
            self.assertEqual(row["dihedral_orbit_count"], 0)
            self.assertEqual(row["histogram_rows"], [])
            self.assertEqual(row["gap_necklaces"], [])
            self.assertTrue(row["is_unsat"])

    def test_proof_interpretation_and_guardrails(self):
        proof = self.payload["exact_proof_certificate"]
        self.assertIn("C(12,5)=792", proof["finite_domain"])
        self.assertIn("12 solutions", proof["center_certificate"])
        self.assertIn("192 solutions", proof["d1_necessity"])
        self.assertIn("36 solutions", proof["d6_necessity"])
        self.assertIn("UNSAT counts [0, 0, 0, 0]", proof["maximality_boundary"])
        self.assertIn("chi_11/shell-label premise", proof["conditional_scope"])

        interpretation = self.payload["interpretation"]
        self.assertIn("locally necessary", interpretation["honest_positive"])
        self.assertIn("does not derive chi_11", interpretation["honest_negative"])
        self.assertIn("Hasse-neighborhood", interpretation["relation_to_previous_probe"])

        ontology = self.payload["ontology_guardrail"]
        self.assertIn("nadsoliton itself", ontology["allowed_reading"])
        self.assertIn("No separate informational layer", ontology["forbidden_reading"])

        hard_limits = "\n".join(self.payload["hard_limits"])
        self.assertIn("No identity K_legacy_ont == K_strict_gate", hard_limits)
        self.assertIn("No theorem derives the chi_11-kernel", hard_limits)
        self.assertIn("conditional finite local necessity certificate", hard_limits)
        self.assertIn("No QW-2191 discharge", hard_limits)
        self.assertIn("No ToE closure", hard_limits)
        self.assertTrue(REPORT_MD.exists())


if __name__ == "__main__":
    unittest.main()
