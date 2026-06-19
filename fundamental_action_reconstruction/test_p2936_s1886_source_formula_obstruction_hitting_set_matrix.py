import json
import subprocess
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
SCRIPT = ROOT / "p2936_s1886_source_formula_obstruction_hitting_set_matrix.py"
OUT = ROOT / "generated" / "p2936_s1886_source_formula_obstruction_hitting_set_matrix.json"
MD = ROOT / "generated" / "p2936_s1886_source_formula_obstruction_hitting_set_matrix.md"


class P2936SourceFormulaObstructionHittingSetMatrixTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        subprocess.run([sys.executable, str(SCRIPT)], cwd=ROOT, check=True, stdout=subprocess.PIPE, text=True)
        cls.payload = json.loads(OUT.read_text(encoding="utf-8"))

    def test_status_and_input(self):
        self.assertEqual(self.payload["status"], "P2936_SOURCE_FORMULA_OBSTRUCTION_HITTING_SET_MATRIX_NO_ACCEPTED_ROUTE")
        self.assertIsNotNone(self.payload["input_hashes"]["P2935"])

    def test_obstruction_certificate(self):
        cert = self.payload["certificate"]
        self.assertEqual(cert["route_count"], 5)
        self.assertEqual(cert["obligation_count"], 5)
        self.assertEqual(cert["accepted_route_count"], 0)
        self.assertEqual(cert["minimal_hitting_set_size"], 1)
        self.assertTrue(cert["minimal_hitting_set_count"] >= 1)
        self.assertTrue(cert["all_current_routes_blocked"])
        self.assertTrue(cert["coordinate_scan_replay_closed"])

    def test_matrix_rows_and_no_closure_flags(self):
        matrix = self.payload["constructed_theoretical_objects"]["obstruction_matrix"]
        self.assertEqual(len(matrix), 5)
        self.assertTrue(all(row["missing_count"] > 0 for row in matrix))
        self.assertTrue(all(not row["accepted_by_P2935_gate"] for row in matrix))
        hits = self.payload["constructed_theoretical_objects"]["minimal_global_hitting_sets_of_missing_obligations"]
        self.assertIn(["strict_nadsoliton_formula_provenance"], hits)
        flags = self.payload["decision"]["negative_export_flags"]
        for key in ["strict_aut_breaking_prime_coordinate_source_law_exported", "strict_prime_log_value_source_exported", "strict_delta_eta_source_exported", "strict_beta_eta_coupling_theorem_exported", "strict_damping_beta_eta_source_packet_exported", "nonproxy_ltotal_exported", "eom_hamiltonian_exported", "bridge_closure_exported", "role_transfer_exported", "toe_closure_exported"]:
            self.assertFalse(flags[key])

    def test_docs_updated(self):
        self.assertIn("P2936/S1886", MD.read_text(encoding="utf-8"))
        self.assertIn("P2936/S1886", (ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md").read_text(encoding="utf-8"))
        self.assertIn("P2936/S1886", (ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md").read_text(encoding="utf-8"))
        self.assertIn("P2936", (REPO / "AGENTS.md").read_text(encoding="utf-8"))


if __name__ == "__main__":
    unittest.main()
