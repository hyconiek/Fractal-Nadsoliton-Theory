import json
import subprocess
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
SCRIPT = ROOT / "p2844_s1794_closure_gate_prime_implicant_obligation_matrix.py"
JSON_PATH = ROOT / "generated" / "p2844_s1794_closure_gate_prime_implicant_obligation_matrix.json"
MD_PATH = ROOT / "generated" / "p2844_s1794_closure_gate_prime_implicant_obligation_matrix.md"


class P2844ClosureGatePrimeImplicantMatrixTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        subprocess.run([sys.executable, str(SCRIPT)], check=True, cwd=ROOT)
        cls.payload = json.loads(JSON_PATH.read_text(encoding="utf-8"))

    def test_status_and_inputs(self):
        self.assertEqual(self.payload["status"], "P2844_CLOSURE_GATE_PRIME_IMPLICANT_OBLIGATION_MATRIX_NO_CLOSURE")
        self.assertEqual(
            self.payload["input_statuses_rechecked"]["P2843"],
            "P2843_PROFESSORIAL_SWOT_CLOSURE_PATH_AUDIT_NO_NEW_CLOSURE",
        )

    def test_current_evidence_atoms_are_narrow(self):
        self.assertEqual(
            set(self.payload["current_evidence_atoms"]),
            {"finite_witness_available", "kernel_split_preserved", "typed_source_domain"},
        )

    def test_target_prime_implicant_matrix(self):
        rows = self.payload["target_prime_implicant_matrix"]
        self.assertEqual(
            set(rows),
            {
                "strict_kernel_completion_bridge",
                "legacy_role_transfer",
                "typed_ltotal_source_coupling",
                "eom_closure",
                "hamiltonian_closure",
                "toe_style_promotion",
            },
        )
        for row in rows.values():
            self.assertFalse(row["currently_closed"])
            self.assertGreater(row["missing_count"], 0)
            self.assertTrue(row["minimal_unlock_sets"])
        self.assertIn("legendre_regular_or_constraint_split", rows["hamiltonian_closure"]["missing_atoms"])
        self.assertIn("nonproxy_ltotal_term", rows["typed_ltotal_source_coupling"]["missing_atoms"])

    def test_candidate_next_move_scores(self):
        scores = {row["candidate_move"]: row for row in self.payload["candidate_next_move_scores"]}
        self.assertIn("unit_bearing_typed_source_coupling_map", scores)
        self.assertTrue(scores["unit_bearing_typed_source_coupling_map"]["admissible_as_next_single_move"])
        self.assertFalse(scores["hamiltonian_legendre_only"]["admissible_as_next_single_move"])
        self.assertEqual(scores["hamiltonian_legendre_only"]["directly_closes_targets"], [])

    def test_acceptance_and_negative_exports(self):
        acceptance = self.payload["acceptance_matrix"]
        self.assertTrue(acceptance["accepted_as_prime_implicant_obligation_matrix"])
        self.assertFalse(acceptance["exports_any_closure"])
        flags = self.payload["decision"]["negative_export_flags"]
        self.assertFalse(flags["bridge_closure_exported"])
        self.assertFalse(flags["role_bearing_ltotal_promoted"])
        self.assertFalse(flags["eom_closure_exported"])
        self.assertFalse(flags["hamiltonian_closure_exported"])
        self.assertFalse(flags["toe_closure_exported"])

    def test_documents_updated(self):
        self.assertIn("P2844/S1794", MD_PATH.read_text(encoding="utf-8"))
        self.assertIn("P2844/S1794", (ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md").read_text(encoding="utf-8"))
        self.assertIn("P2844/S1794", (ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md").read_text(encoding="utf-8"))
        self.assertIn("P2844", (REPO / "AGENTS.md").read_text(encoding="utf-8"))


if __name__ == "__main__":
    unittest.main()
