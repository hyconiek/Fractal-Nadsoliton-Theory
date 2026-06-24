import json, subprocess, sys, unittest
from pathlib import Path
ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
SCRIPT = ROOT / "p3077_s2027_second_order_lift_obstruction_table.py"
OUT = ROOT / "generated" / "p3077_s2027_second_order_lift_obstruction_table.json"
MD = ROOT / "generated" / "p3077_s2027_second_order_lift_obstruction_table.md"

class P3077SecondOrderLiftObstructionTableTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        subprocess.run([sys.executable, str(SCRIPT)], cwd=ROOT, check=True, stdout=subprocess.PIPE, text=True)
        cls.payload = json.loads(OUT.read_text(encoding="utf-8"))

    def test_status_and_input_hash(self):
        self.assertEqual(self.payload["status"], "P3077_FORMAL_HAMILTONIAN_LIFT_EXISTS_INTERNAL_SECOND_ORDER_SOURCE_OBSTRUCTED")
        self.assertIsNotNone(self.payload["input_hashes"]["P3076"])

    def test_matrix_counts(self):
        cert = self.payload["finite_certificate"]
        self.assertEqual(cert["content_grep_lanes"], 4)
        self.assertGreater(cert["content_grep_hits"], 0)
        self.assertEqual(cert["p3076_z12_modes"], 12)
        self.assertEqual(cert["p3076_accepted_lightlike_branch_rows"], 0)
        self.assertEqual(cert["candidate_lifts"], 4)
        self.assertEqual(cert["premise_gates"], 6)
        self.assertEqual(cert["premise_gate_rows"], 24)
        self.assertEqual(cert["mode_lift_rows"], 48)
        self.assertEqual(cert["formal_wave_compatible_rows"], 11)
        self.assertEqual(cert["formal_hamiltonian_wave_rows"], 11)
        self.assertEqual(cert["accepted_internal_second_order_wave_rows"], 0)
        self.assertEqual(cert["internally_sourced_required_wave_gate_rows"], 0)
        self.assertEqual(cert["satisfied_proof_obligations"], 5)

    def test_aggregate_and_formal_rows(self):
        rows = self.payload["constructed_theoretical_objects"]["candidate_aggregate_certificate"]
        by_id = {row["candidate_lift"]: row for row in rows}
        self.assertEqual(by_id["formal_hamiltonian_dirichlet_lift"]["formal_wave_compatible_rows"], 11)
        self.assertEqual(by_id["formal_hamiltonian_dirichlet_lift"]["accepted_internal_second_order_wave_rows"], 0)
        self.assertEqual(by_id["gradient_flow_only_baseline"]["formal_wave_compatible_rows"], 0)
        reps = self.payload["constructed_theoretical_objects"]["representative_formal_hamiltonian_rows"]
        self.assertTrue(reps)
        self.assertTrue(all(row["formal_wave_compatible_row"] for row in reps))
        self.assertTrue(all(not row["accepted_internal_second_order_wave_row"] for row in reps))

    def test_flags_recommendation_and_docs(self):
        for value in self.payload["decision"]["negative_export_flags"].values():
            self.assertFalse(value)
        self.assertTrue(self.payload["decision"]["positive_scoped_flags"]["formal_hamiltonian_dirichlet_lift_constructed"])
        self.assertIn("momentum/symplectic-source audit", self.payload["decision"]["next_honest_step"])
        self.assertIn("P3077/S2027", MD.read_text(encoding="utf-8"))
        self.assertIn("P3077/S2027", (ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md").read_text(encoding="utf-8"))
        self.assertIn("P3077/S2027", (ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md").read_text(encoding="utf-8"))
        self.assertIn("P3077", (REPO / "AGENTS.md").read_text(encoding="utf-8"))

if __name__ == "__main__":
    unittest.main()
